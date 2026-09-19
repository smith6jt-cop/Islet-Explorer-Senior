#!/usr/bin/env python3
"""End-to-end smoke test for the TissUUmaps viewer path.

Builds a `.tmap` project, boots a TissUUmaps server against it, asserts every
route the app depends on, and (optionally) drives a headless browser to confirm
the layers, islet regions and single-cell markers actually render.

Run it two ways:

    # Synthetic slides — no cohort data needed. Reproduces the batch-1 (59-ch,
    # SST@35) and batch-2 (58-ch, SST@1) panel layouts, so the by-name channel
    # resolution is genuinely exercised.
    python scripts/senior/tissuumaps_smoke_test.py

    # A real donor slide — this is the one that matters. Prints the true
    # conversion time and disk ratio for the channel pyramids.
    python scripts/senior/tissuumaps_smoke_test.py \
        --image /data/islet_ome_tiff/6539.ome.tiff \
        --regions data/Islets \
        --markers data/tissuumaps_inputs/markers/6539.csv

Requires: tissuumaps, pyvips, tifffile (+ playwright for --browser).
"""

from __future__ import annotations

import argparse
import gzip
import json
import math
import os
import random
import re
import socket
import subprocess
import sys
import time
import urllib.request
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
BUILDER = REPO / "scripts" / "senior" / "build_tissuumaps_project.py"
CHANNEL_NAMES = REPO / "apps" / "islet" / "Channel_names"

PASS, FAIL = "PASS", "FAIL"
results: list[tuple[str, str, str]] = []


def check(label: str, ok: bool, detail: str = "") -> bool:
    results.append((label, PASS if ok else FAIL, detail))
    print(f"  {PASS if ok else FAIL}  {label}" + (f"  [{detail}]" if detail else ""), flush=True)
    return ok


def log(msg: str) -> None:
    print(f"\n=== {msg}", flush=True)


# ---------------------------------------------------------------------------
# Synthetic slides
# ---------------------------------------------------------------------------


def panel_names(order: str) -> list[str]:
    names = [ln.split(" (C")[0].strip() for ln in CHANNEL_NAMES.read_text().splitlines()
             if ln.strip()]
    if order == "batch2":
        # Batch-2 slides drop CD38 (58 channels) and put SST at index 1.
        names = [n for n in names if n != "CD38"]
        names.remove("SST")
        names.insert(1, "SST")
    return names


def make_slide(path: Path, order: str, size: int = 1024) -> list[str]:
    import numpy as np
    import tifffile

    names = panel_names(order)
    rng = np.random.default_rng(0)
    yy, xx = np.mgrid[0:size, 0:size]
    data = np.zeros((len(names), size, size), np.uint16)
    for i in range(len(names)):
        cx, cy = 300 + (i * 37) % 400, 300 + (i * 53) % 400
        r = 60 + (i % 5) * 20
        blob = np.exp(-(((xx - cx) ** 2 + (yy - cy) ** 2) / (2.0 * r * r)))
        # Background is deliberately SKEWED, not just noisy: a symmetric normal
        # floor cannot reproduce the failure the stretch estimator has to
        # survive. Real autofluorescence puts a long upper tail above the
        # background peak, and it was that tail -- counted as signal -- which
        # dragged the white point below the marker and blew every rare channel
        # out. With a normal-only floor the old and new estimators score the
        # same and the regression check below is worthless.
        floor = rng.normal(200, 60, (size, size)) + rng.gamma(1.2, 70, (size, size))
        # Every third channel is a RARE marker: bright, but over few pixels, so
        # its population sits under the background tail's. That is the case that
        # broke, and it has to be in the fixture.
        if i % 3 == 0:
            mask = blob > 0.75
            amp = np.zeros((size, size))
            amp[mask] = 12000 + 300 * i
        else:
            amp = blob * (8000 + 400 * i)
        data[i] = np.clip(amp + floor, 0, 65535)

    path.parent.mkdir(parents=True, exist_ok=True)
    with tifffile.TiffWriter(str(path), bigtiff=True) as tif:
        opts = dict(tile=(256, 256), photometric="minisblack")
        tif.write(data, subifds=2, metadata={"axes": "CYX", "Channel": {"Name": names}}, **opts)
        for lvl in (2, 4):  # SubIFD pyramid levels, as raw2ometiff writes them
            tif.write(data[:, ::lvl, ::lvl], subfiletype=1, **opts)
    return names


def make_structures(path: Path, n: int = 300) -> int:
    """Vessel/nerve/lymphatic GeoJSON, in the same resolution-#1 pixel space.

    Small, many, and multi-class -- the shape of the real exports rather than
    their size (those run to ~200k features per donor). The point of the fixture
    is to prove the SECOND region layer is emitted, titled and left un-auto-loaded;
    proving 200k features render acceptably needs a real donor, not a synthetic one.
    """
    random.seed(7)
    classes = ["Capillary", "Nerve", "Lymphatic"]
    feats = []
    for k in range(n):
        cx, cy, r = random.uniform(20, 460), random.uniform(20, 400), random.uniform(3, 9)
        ring = [[cx + r * math.cos(t * math.pi / 4), cy + r * math.sin(t * math.pi / 4)]
                for t in range(8)]
        ring.append(ring[0])
        cls = classes[k % len(classes)]
        feats.append({
            "type": "Feature", "id": f"{cls}_{k}",
            "geometry": {"type": "Polygon", "coordinates": [ring]},
            "properties": {"objectType": "annotation", "classification": {"name": cls}},
        })
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps({"type": "FeatureCollection", "features": feats}))
    return len(feats)


def make_aux(geojson: Path, csv: Path, um_per_px_res1: float = 0.5078) -> tuple[int, int]:
    """QuPath-style islet GeoJSON (resolution-#1 px) + per-cell CSV (µm)."""
    random.seed(0)
    feats, rows = [], ["X_centroid,Y_centroid,phenotype,cell_region,islet_name"]
    phenos = ["Beta cell", "Alpha cell", "Delta cell", "Macrophage", "CD8a Tcell", "Endothelial"]
    for k in range(6):
        cx, cy, r = 80 + k * 70, 120 + (k % 3) * 130, 22 + k * 3
        ring = [[cx + r * math.cos(t * math.pi / 12), cy + r * math.sin(t * math.pi / 12)]
                for t in range(24)]
        ring.append(ring[0])
        feats.append({
            "type": "Feature", "id": f"Islet_{k+1}",
            "geometry": {"type": "Polygon", "coordinates": [ring]},
            "properties": {"objectType": "annotation", "name": f"Islet_{k+1}",
                           "classification": {"name": "Islet"}},
        })
        for _ in range(400):
            a, d = random.random() * 2 * math.pi, r * math.sqrt(random.random()) * 1.4
            rows.append(",".join([
                f"{(cx + d * math.cos(a)) * um_per_px_res1:.2f}",
                f"{(cy + d * math.sin(a)) * um_per_px_res1:.2f}",
                random.choice(phenos), "core" if d < r else "peri", f"Islet_{k+1}",
            ]))

    geojson.parent.mkdir(parents=True, exist_ok=True)
    blob = json.dumps({"type": "FeatureCollection", "features": feats})
    if geojson.suffix == ".gz":
        with gzip.open(geojson, "wt") as fh:
            fh.write(blob)
    else:
        geojson.write_text(blob)
    csv.write_text("\n".join(rows) + "\n")
    return len(feats), len(rows) - 1


def case_of(image: Path) -> str:
    """Case id the builder will derive from a slide filename."""
    return re.sub(r"\.ome$", "", image.stem, flags=re.I)


def image_size(image: Path) -> tuple[int, int]:
    import tifffile

    with tifffile.TiffFile(str(image)) as tif:
        return int(tif.pages[0].imagewidth), int(tif.pages[0].imagelength)


def find_region_source(spec: Path | None, case_id: str) -> Path | None:
    """The GeoJSON the builder would have read, given a file or a directory.

    Covers the same extensions the builder does. A helper that only knew about
    plain ``.geojson`` would quietly skip the scale assertion below on a gzipped
    or ``.json`` export -- i.e. it would go green precisely where the coordinate
    handling is least exercised.
    """
    if spec is None or not spec.exists():
        return None
    if spec.is_file():
        return spec
    bounded = re.compile(rf"(?<![0-9]){re.escape(case_id)}(?![0-9])")
    hits = sorted({f for pat in ("*.geojson", "*.geojson.gz", "*.json", "*.json.gz")
                   for f in spec.glob(pat) if bounded.search(f.name)})
    return hits[0] if hits else None


def pyramid_level_means(tif: Path) -> list[float]:
    """Mean pixel value of each pyramid level, read the way the server reads it.

    Via libvips (libtiff), i.e. the same decoder family OpenSlide uses to cut
    tiles -- so this measures what viewers actually get, not what a second
    decoder thinks the file says.

    Uses vips' own ``avg()``, which streams: pulling each level through
    ``write_to_memory()`` instead would allocate the whole raster in Python,
    and level 0 of a real donor channel is tens of thousands of pixels per side
    -- hundreds of MB per allocation, in the run this check exists to protect.
    ``avg()`` also averages over all bands and any pixel format, so it carries
    no assumption that the file is single-band uint8.
    """
    import pyvips

    means, page = [], 0
    while True:
        try:
            v = pyvips.Image.new_from_file(str(tif), page=page)
        except Exception:  # noqa: BLE001 - ran off the end of the directory list
            break
        means.append(float(v.avg()))
        if v.width <= 256 and v.height <= 256:
            break
        page += 1
    return means


def read_geojson(path: Path):
    """Load a GeoJSON, gzipped or not."""
    if path.suffix == ".gz":
        with gzip.open(path, "rt") as fh:
            return json.load(fh)
    return json.loads(path.read_text())


def geojson_max(gj) -> tuple[float, float] | None:
    """Largest x and y over every coordinate. Deliberately re-implemented here
    rather than imported from the builder, so the test can catch a bug in the
    builder's own bounds walk."""
    best = [float("-inf"), float("-inf")]

    def walk(c):
        if c and isinstance(c[0], (int, float)):
            best[0], best[1] = max(best[0], c[0]), max(best[1], c[1])
        else:
            for sub in c:
                walk(sub)

    feats = gj.get("features", [gj]) if isinstance(gj, dict) else gj
    for feat in feats:
        geom = (feat or {}).get("geometry")
        if geom and geom.get("coordinates") is not None:
            walk(geom["coordinates"])
    return None if best[0] == float("-inf") else (best[0], best[1])


# ---------------------------------------------------------------------------
# Server
# ---------------------------------------------------------------------------


def free_port() -> int:
    with socket.socket() as s:
        s.bind(("127.0.0.1", 0))
        return s.getsockname()[1]


def wait_for(url: str, timeout: float = 60.0) -> bool:
    deadline = time.time() + timeout
    while time.time() < deadline:
        try:
            with urllib.request.urlopen(url, timeout=2) as r:
                if r.status == 200:
                    return True
        except Exception:
            time.sleep(0.5)
    return False


def http_status(url: str) -> tuple[int, int]:
    """(status, bytes). Returns (0, 0) on connection failure."""
    try:
        with urllib.request.urlopen(url, timeout=30) as r:
            return r.status, len(r.read())
    except urllib.error.HTTPError as e:
        return e.code, 0
    except Exception:
        return 0, 0


def browser_check(url: str) -> dict | None:
    try:
        from playwright.sync_api import sync_playwright
    except ImportError:
        print("  SKIP browser check (playwright not installed)")
        return None

    exe = os.environ.get("CHROMIUM_PATH", "/opt/pw-browsers/chromium")
    launch = {"args": ["--no-sandbox", "--use-gl=swiftshader", "--enable-unsafe-swiftshader"]}
    if Path(exe).exists():
        launch["executable_path"] = exe

    with sync_playwright() as p:
        b = p.chromium.launch(**launch)
        pg = b.new_page(viewport={"width": 1500, "height": 950})
        pg.goto(url, wait_until="networkidle", timeout=120000)
        pg.wait_for_timeout(10000)
        info = pg.evaluate("""() => {
            const vis = {};
            (tmapp.layers || []).forEach((l, i) => {
                const cb = document.getElementById('visible-layer-' + i);
                if (cb && cb.checked) vis[l.name] = i;
            });
            return {
                layers: (tmapp.layers || []).length,
                worldItems: tmapp['ISS_viewer'] ? tmapp['ISS_viewer'].world.getItemCount() : -1,
                visible: vis,
                composite: filterUtils._compositeMode,
                regions: Object.keys(regionUtils._regions || {}).length,
                // Only VISIBLE modals matter, and only the blocking ones: the
                // transient "Loading image, please wait..." notification lives
                // in the same container. The one that breaks region autoLoad is
                // the layer-index prompt (suppressed by properties.collectionIndex).
                modals: [...document.querySelectorAll('.modal-body')]
                    .filter(e => e.offsetParent !== null)
                    .map(e => e.innerText.slice(0, 60))
                    .filter(t => !/loading/i.test(t)),
            };
        }""")
        btns = pg.locator("#divMarkersDownloadButtons button")
        if btns.count():
            btns.nth(0).click()
            pg.wait_for_timeout(9000)
        info["markerGroups"] = pg.evaluate("""() => {
            const d = Object.values(dataUtils.data || {})[0];
            return d ? Object.keys(d._groupgarden || {}).length : 0;
        }""")
        shot = Path(os.environ.get("SMOKE_SHOT", "")) if os.environ.get("SMOKE_SHOT") else None
        if shot:
            pg.screenshot(path=str(shot))
            info["screenshot"] = str(shot)
        b.close()
    return info


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--image", type=Path, default=None,
                    help="Real donor OME-TIFF (default: synthesize batch-1 + batch-2 slides)")
    ap.add_argument("--regions", type=Path, default=None, help="QuPath GeoJSON (real runs)")
    ap.add_argument("--structures", type=Path, default=None,
                    help="Vessel/nerve/lymphatic GeoJSON (real runs)")
    ap.add_argument("--markers", type=Path, default=None, help="Per-cell CSV (real runs)")
    ap.add_argument("--work-dir", type=Path, default=Path("/tmp/tissuumaps_smoke"),
                    help="Scratch directory for slides + project output")
    ap.add_argument("--channels", default=None,
                    help="Limit extraction to these markers (real runs — much faster)")
    ap.add_argument("--no-browser", action="store_true", help="Skip the headless-browser check")
    ap.add_argument("--keep", action="store_true", help="Leave the server running at the end")
    args = ap.parse_args()

    work = args.work_dir.resolve()
    projects = work / "projects"
    work.mkdir(parents=True, exist_ok=True)

    # -- slides -------------------------------------------------------------
    if args.image:
        slides = [(args.image.resolve(), args.regions, args.markers, args.structures)]
        log(f"Using real slide {args.image}")
    else:
        log("Building synthetic slides (batch-1 59ch SST@35, batch-2 58ch SST@1)")
        # Named, located and compressed the way the real exports are (a
        # directory of per-donor files, often gzipped), so --regions exercises
        # the directory lookup and the gzip path rather than a plain file the
        # builder rarely sees in practice.
        gj = work / "aux" / "results__6539_IO60panc2nd_resolution__1.geojson.gz"
        csv = work / "aux" / "cells.csv"
        n_islets, n_cells = make_aux(gj, csv)
        # Structures live in their OWN directory, as they do in reality (islet
        # exports in data/Islets, structures written by build_viewer_assets.py
        # into app_data/json). Co-locating them would make the case-id lookup
        # ambiguous -- `6539_structures.geojson` also matches a search for 6539 --
        # and the builder would fail closed and attach neither.
        struct = work / "structures" / "6539_structures.geojson"
        n_struct = make_structures(struct)
        print(f"  aux: {n_islets} islets, {n_struct} structures, {n_cells} cells")
        slides = []
        for case, order in (("6539", "batch1"), ("6450", "batch2")):
            path = work / "slides" / f"{case}.ome.tiff"
            if not path.exists():
                names = make_slide(path, order)
                print(f"  {case}: {len(names)} channels, SST@{names.index('SST')}")
            slides.append((path, gj.parent if case == "6539" else None,
                           csv if case == "6539" else None,
                           struct.parent if case == "6539" else None))

    # -- build projects -----------------------------------------------------
    log("Building .tmap projects")
    for image, regions, markers, structures in slides:
        cmd = [sys.executable, str(BUILDER), "--image", str(image),
               "--out-dir", str(projects)]
        if regions:
            cmd += ["--regions", str(regions)]
        if structures:
            cmd += ["--structures", str(structures)]
        if markers:
            cmd += ["--markers", str(markers)]
        if args.channels:
            cmd += ["--channels", args.channels]
        t0 = time.time()
        proc = subprocess.run(cmd, capture_output=True, text=True)
        tail = [ln for ln in proc.stdout.splitlines()
                if any(k in ln for k in ("channels:", "regions:", "markers:", "default-on",
                                         "schema", "-> index", "note:", "WARNING",
                                         "ERROR"))]
        for ln in tail:
            print("   ", ln)
        if not check(f"build {image.name}", proc.returncode == 0,
                     f"{time.time() - t0:.0f}s"):
            print(proc.stderr[-2000:])
            return 1

    cases = sorted(p.stem for p in projects.glob("*.tmap"))
    check("projects written", bool(cases), ", ".join(cases))

    primary = cases[0] if len(cases) == 1 else "6539"
    tmap = json.loads((projects / f"{primary}.tmap").read_text())
    names = [layer["name"] for layer in tmap["layers"]]
    visible = [tmap["layers"][int(i)]["name"] for i, v in tmap["layerVisibilities"].items() if v]
    check("layers named after markers", not any(n.startswith("Channel_") for n in names),
          f"{len(names)} layers, e.g. {names[:3]}")
    check("default-on channels resolved by name",
          set(visible) <= {"DAPI", "INS", "SST", "GCG"} and len(visible) > 0,
          ", ".join(sorted(visible)))
    # "lighter" is canvas' additive operator, and additive-with-clamp is exactly
    # how QuPath composites fluorescence channels -- verified in Chromium against
    # a clamped-sum reference: of source-over/lighter/screen/multiply/darken/
    # overlay/color-dodge/difference/luminosity, only lighter reproduced it.
    # TissUUmaps' own default is source-over, which shows the top layer ONLY.
    check("additive blending (QuPath-equivalent)", tmap["compositeMode"] == "lighter")

    # Every layer must come up in its own visible colour. White is as bad as
    # black here: under additive blending it lands on all three components, so
    # switching on a fifth marker both washes out the others and is
    # indistinguishable from switching on any other.
    colors = {layer["name"]: tmap["layerFilters"][str(i)][0]["value"]
              for i, layer in enumerate(tmap["layers"])}
    unusable = {n: c for n, c in colors.items() if c in ("100,100,100", "0,0,0")}
    check("no layer starts white or black", not unusable, str(unusable) or f"{len(colors)} layers")
    lit = all(max(int(v) for v in c.split(",")) >= 50 for c in colors.values())
    check("every layer colour is fully lit", lit)

    # Reduced pyramid levels must stay photometrically consistent with level 0.
    # Downsampling averages, so a level's mean can drift a little (edges, and
    # clipping at the white point) but cannot climb: libvips' default horizontal
    # predictor silently corrupted every reduced level, which reads as the
    # background lifting off black as you zoom out while full zoom looks
    # perfect. Nothing in the .tmap or the HTTP routes can catch that -- only
    # the pixels can.
    for ch_path in sorted((projects / "channels" / primary).glob("*.tif"))[:4]:
        means = pyramid_level_means(ch_path)
        if len(means) < 2:
            check(f"{ch_path.name} has reduced levels", False, "single level")
            continue
        base = max(1.0, means[0])
        worst = max(abs(m - means[0]) / base for m in means[1:])
        check(f"{ch_path.name} levels consistent with level 0",
              worst < 0.35, "means " + ", ".join(f"{m:.1f}" for m in means))

    # Regions are exported from QuPath at a downsampled pyramid level, so the
    # builder rescales them. Getting that wrong is a silent misalignment, not an
    # error, which is exactly what a test is for: the outlines must land inside
    # the full-resolution image and must have actually moved.
    # The stretch must not blow the markers out, and this has to be shown BOTH
    # ways: a check that only asserts the good case would also pass if the
    # fixture were too easy to fail on, which is exactly what happened before
    # the synthetic background grew an asymmetric tail.
    hit = next((s for s in slides if case_of(s[0]) == primary), None)
    if hit is not None:
        def report(floor_k: str) -> tuple[str, str]:
            """(verdict, detail) from a --report run, or ('error', why).

            The builder's exit status and stderr are surfaced rather than
            dropped: a report that failed to run produces no verdict line, and
            without this the assertions below would report a stretch problem
            when the real fault was a missing dependency or a CLI regression.
            """
            cmd = [sys.executable, str(BUILDER), "--image", str(hit[0]),
                   "--out-dir", str(projects), "--report", "--signal-floor-k", floor_k]
            if args.channels:
                cmd += ["--channels", args.channels]
            proc = subprocess.run(cmd, capture_output=True, text=True)
            line = next((ln for ln in proc.stdout.splitlines() if "VERDICT:" in ln), None)
            if proc.returncode != 0 or line is None:
                return "error", (proc.stderr.strip().splitlines() or ["no output"])[-1][:160]
            verdict, _, detail = line.split("VERDICT:", 1)[1].strip().partition(" - ")
            return verdict, detail

        # Marker coordinates are trusted to be µm and multiplied by 1/µm-per-px.
        # If they were already pixels every cell lands ~4x out -- which looks
        # unremarkable zoomed out and is only obvious at cell-level zoom. Prove
        # the builder refuses rather than attaching them, by feeding it the same
        # CSV with the coordinates pre-multiplied so they overflow the slide.
        if hit[2]:
            import csv as _csv

            # Into its OWN out-dir, and one channel: this run writes a .tmap for
            # the same case id, so sharing the project directory would overwrite
            # the real one -- without its regions, because this invocation does
            # not pass them.
            probe_dir = projects.parent / "marker_probe"
            bad = probe_dir / "_bad_markers.csv"
            bad.parent.mkdir(parents=True, exist_ok=True)
            with open(hit[2], newline="") as fh, open(bad, "w", newline="") as out:
                r = _csv.reader(fh)
                w = _csv.writer(out)
                head = next(r)
                w.writerow(head)
                for row in r:
                    row[0], row[1] = str(float(row[0]) * 40), str(float(row[1]) * 40)
                    w.writerow(row)
            proc = subprocess.run(
                [sys.executable, str(BUILDER), "--image", str(hit[0]),
                 "--out-dir", str(probe_dir), "--markers", str(bad), "--channels", "DAPI"],
                capture_output=True, text=True)
            check("out-of-bounds markers are refused, not attached",
                  proc.returncode == 3 and "NOT attaching" in proc.stdout,
                  f"exit {proc.returncode}")

        good_v, good_d = report("20")
        old_v, old_d = report("3")
        check("stretch does not blow out the markers", good_v == "ok",
              good_d if good_v != "error" else f"report failed: {good_d}")
        # Without this second direction the first check would also pass on a
        # fixture too easy to fail on -- which is what the synthetic slides were
        # until they grew an asymmetric background tail and rare bright markers.
        check("that check can fail (old floor is caught)", old_v == "blown-out",
              old_d if old_v != "error" else f"report failed: {old_d}")
    # Islets must auto-load and structures must not: the vessel/nerve/lymphatic
    # exports run to ~200k features per donor, and auto-loading them would make
    # every donor pay that cost just to look at an islet.
    region_layers = {r["title"]: r for r in tmap.get("regionFiles", [])}
    auto = [t for t, r in region_layers.items() if r.get("autoLoad")]
    if region_layers:
        check("islets auto-load, nothing else does",
              auto == [t for t in region_layers if t.startswith("Islets")],
              f"{len(region_layers)} layer(s); auto-loading: {auto or 'none'}")

    # Without per-class colours every structure class draws the same, so
    # vessels, nerves and lymphatics overlaid together are indistinguishable --
    # which is the whole reason for showing them.
    struct_layer = next((r for t, r in region_layers.items() if t.startswith("Structures")), None)
    if struct_layer:
        cls_colors = {}
        for f in read_geojson(projects / struct_layer["path"]).get("features", []):
            cls = (f.get("properties") or {}).get("classification") or {}
            if cls.get("name") and cls.get("color"):
                cls_colors[cls["name"]] = tuple(cls["color"])
        check("structure classes get distinct colours",
              len(cls_colors) > 1 and len(set(cls_colors.values())) == len(cls_colors),
              ", ".join(f"{k}={v}" for k, v in sorted(cls_colors.items())) or "none coloured")

    for rel in (r["path"] for r in region_layers.values()):
        if hit is None:
            break
        img_w, img_h = image_size(hit[0])
        out_max = geojson_max(read_geojson(projects / rel))
        check(f"{rel} lies inside the {img_w}x{img_h} image",
              out_max is not None and out_max[0] <= img_w and out_max[1] <= img_h,
              f"max {out_max[0]:.0f},{out_max[1]:.0f}" if out_max else "no geometry")
        # Compare each written layer against the source it actually came from --
        # islets to hit[1], structures to hit[3]. Measuring both against the
        # islet export would silently pass whatever the structures did.
        stem = Path(rel).stem
        spec, want = (hit[3], stem) if stem.endswith("_structures") else (hit[1], primary)
        src_file = find_region_source(spec, want)
        src_max = geojson_max(read_geojson(src_file)) if src_file else None
        if src_max and src_max[0] and out_max:
            factor = out_max[0] / src_max[0]
            check(f"{stem}: region scale is a clean pyramid factor",
                  any(abs(factor - c) < 0.01 for c in (0.25, 0.5, 1, 2, 4, 8)), f"x{factor:.2f}")

    if len(cases) > 1:  # synthetic run: prove the panels differ but names resolve
        other = json.loads((projects / "6450.tmap").read_text())
        idx = {layer["name"]: i for i, layer in enumerate(other["layers"])}
        check("batch-2 panel resolves independently",
              idx.get("SST") != names.index("SST") if "SST" in idx else False,
              f"batch-1 SST@{names.index('SST')} vs batch-2 SST@{idx.get('SST')}")

        # The reason layer colours are hashed from the marker NAME rather than
        # its index: the two panels order markers differently, so an index-based
        # colour would render the same marker differently on different donors.
        # This is the one check that can prove it, because it is the only place
        # two panels exist side by side.
        other_colors = {layer["name"]: other["layerFilters"][str(i)][0]["value"]
                        for i, layer in enumerate(other["layers"])}
        shared = set(colors) & set(other_colors)
        drifted = {n for n in shared if colors[n] != other_colors[n]}
        check("colours are stable across panels", not drifted,
              f"{len(shared)} shared markers"
              + (f", drifted: {sorted(drifted)}" if drifted else ""))

    # -- server -------------------------------------------------------------
    # Report which build is under test. The deployment runs our fork
    # (smith6jt-cop/TissUUmaps), so a run against a pip-installed upstream is
    # testing something other than what ships -- make that visible in the log
    # rather than leaving it to be inferred.
    # This banner is diagnostics, never a gate: importing tissuumaps drags in
    # pyvips/h5py/openslide and can fail with things that are not ImportError
    # (a mismatched libhdf5 raises ValueError, for one -- see § 6 of
    # docs/tissuumaps_evaluation.md). The server runs out of process, so a
    # failure here says nothing about whether the run can proceed.
    try:
        import tissuumaps as _tm

        # There is no __version__; the package exposes tissuumaps_version and
        # stashes the same string in app.config["VERSION"].
        version = getattr(_tm, "tissuumaps_version", None)
        if not version:
            version = getattr(_tm, "app", None) and _tm.app.config.get("VERSION")
        log(f"TissUUmaps {version or '?'} from {Path(_tm.__file__).parent}")
    except Exception as exc:  # noqa: BLE001 - banner must never fail the run
        log(f"could not identify the installed TissUUmaps ({type(exc).__name__}: {exc}); "
            f"using tissuumaps_server from PATH")

    port = free_port()
    base = f"http://127.0.0.1:{port}"
    # Prefer the entrypoint that sits beside the interpreter running this test:
    # running it from a venv without activating the venv is the normal case, and
    # falling through to PATH there would either fail outright or -- worse --
    # silently test a different TissUUmaps than the banner above just reported.
    exe = Path(sys.executable).parent / "tissuumaps_server"
    server_cmd = str(exe) if exe.exists() else "tissuumaps_server"
    log(f"Starting {server_cmd} on {port}")
    server = subprocess.Popen(
        [server_cmd, "-l", "127.0.0.1", "-p", str(port), str(projects)],
        stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT,
    )
    try:
        if not check("server responds to /ping", wait_for(f"{base}/ping")):
            print("  (if this fails on import, see docs/tissuumaps_evaluation.md § 6)")
            return 1

        layer0 = tmap["layers"][0]["tileSource"]
        routes = [
            (f"/{primary}.tmap", 200),
            ("/" + layer0, 200),
            ("/" + layer0.replace(".dzi", "_files/10/0_0.jpeg"), 200),
            ("/nonexistent.tmap", 404),
        ]
        for rel in (r["path"] for r in tmap.get("regionFiles", [])):
            routes.append((f"/{rel}", 200))
        for rel in (m["path"] for m in tmap.get("markerFiles", [])):
            routes.append((f"/{rel}", 200))

        for rel, want in routes:
            status, size = http_status(base + rel)
            check(f"GET {rel}", status == want, f"{status}, {size}B")

        # -- browser --------------------------------------------------------
        if not args.no_browser:
            log("Headless browser check")
            info = browser_check(f"{base}/{primary}.tmap")
            if info:
                check("all layers mounted in the viewer",
                      info["worldItems"] == info["layers"], f"{info['worldItems']} items")
                check("only the marker channels visible",
                      set(info["visible"]) <= {"DAPI", "INS", "SST", "GCG"},
                      ", ".join(info["visible"]))
                check("no blocking modal", not info["modals"],
                      "; ".join(info["modals"])[:80])
                if tmap.get("regionFiles"):
                    check("islet regions auto-loaded", info["regions"] > 0,
                          f"{info['regions']} regions")
                if tmap.get("markerFiles"):
                    check("cell markers load and group by phenotype",
                          info["markerGroups"] > 0, f"{info['markerGroups']} groups")

        if args.keep:
            print(f"\nServer left running at {base}/{primary}.tmap (Ctrl-C to stop)")
            server.wait()
    finally:
        if not args.keep:
            server.terminate()
            try:
                server.wait(timeout=10)
            except subprocess.TimeoutExpired:
                server.kill()

    failed = [r for r in results if r[1] == FAIL]
    log(f"{len(results) - len(failed)}/{len(results)} checks passed")
    for label, _, detail in failed:
        print(f"  FAILED: {label} {detail}")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
