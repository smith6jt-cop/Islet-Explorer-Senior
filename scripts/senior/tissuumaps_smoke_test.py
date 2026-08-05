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
        --regions data/app_data/json/6539.geojson \
        --markers data/app_data/cells/6539_tissue.csv

Requires: tissuumaps, pyvips, tifffile (+ playwright for --browser).
"""

from __future__ import annotations

import argparse
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
CHANNEL_NAMES = REPO / "app" / "shiny_app" / "Channel_names"

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
        data[i] = np.clip(blob * (8000 + 400 * i) + rng.normal(200, 60, (size, size)), 0, 65535)

    path.parent.mkdir(parents=True, exist_ok=True)
    with tifffile.TiffWriter(str(path), bigtiff=True) as tif:
        opts = dict(tile=(256, 256), photometric="minisblack")
        tif.write(data, subifds=2, metadata={"axes": "CYX", "Channel": {"Name": names}}, **opts)
        for lvl in (2, 4):  # SubIFD pyramid levels, as raw2ometiff writes them
            tif.write(data[:, ::lvl, ::lvl], subfiletype=1, **opts)
    return names


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
    geojson.write_text(json.dumps({"type": "FeatureCollection", "features": feats}))
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
    """The GeoJSON the builder would have read, given a file or a directory."""
    if spec is None or not spec.exists():
        return None
    if spec.is_file():
        return spec
    return next(iter(sorted(spec.glob(f"*{case_id}*.geojson"))), None)


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
        slides = [(args.image.resolve(), args.regions, args.markers)]
        log(f"Using real slide {args.image}")
    else:
        log("Building synthetic slides (batch-1 59ch SST@35, batch-2 58ch SST@1)")
        # Named and located the way the real exports are (a directory of
        # per-donor files), so --regions exercises the directory lookup rather
        # than a path the builder never sees in practice.
        gj = work / "aux" / "results__6539_IO60panc2nd_resolution__1.geojson"
        csv = work / "aux" / "cells.csv"
        n_islets, n_cells = make_aux(gj, csv)
        print(f"  aux: {n_islets} islets, {n_cells} cells")
        slides = []
        for case, order in (("6539", "batch1"), ("6450", "batch2")):
            path = work / "slides" / f"{case}.ome.tiff"
            if not path.exists():
                names = make_slide(path, order)
                print(f"  {case}: {len(names)} channels, SST@{names.index('SST')}")
            slides.append((path, gj.parent if case == "6539" else None,
                           csv if case == "6539" else None))

    # -- build projects -----------------------------------------------------
    log("Building .tmap projects")
    for image, regions, markers in slides:
        cmd = [sys.executable, str(BUILDER), "--image", str(image),
               "--out-dir", str(projects)]
        if regions:
            cmd += ["--regions", str(regions)]
        if markers:
            cmd += ["--markers", str(markers)]
        if args.channels:
            cmd += ["--channels", args.channels]
        t0 = time.time()
        proc = subprocess.run(cmd, capture_output=True, text=True)
        tail = [ln for ln in proc.stdout.splitlines()
                if any(k in ln for k in ("channels:", "regions:", "markers:", "default-on",
                                         "schema", "-> index", "WARNING", "ERROR"))]
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
    check("additive blending", tmap["compositeMode"] == "lighter")

    # Regions are exported from QuPath at a downsampled pyramid level, so the
    # builder rescales them. Getting that wrong is a silent misalignment, not an
    # error, which is exactly what a test is for: the outlines must land inside
    # the full-resolution image and must have actually moved.
    hit = next((s for s in slides if case_of(s[0]) == primary), None)
    for rel in (r["path"] for r in tmap.get("regionFiles", [])):
        if hit is None:
            break
        img_w, img_h = image_size(hit[0])
        out_max = geojson_max(json.loads((projects / rel).read_text()))
        check(f"{rel} lies inside the {img_w}x{img_h} image",
              out_max is not None and out_max[0] <= img_w and out_max[1] <= img_h,
              f"max {out_max[0]:.0f},{out_max[1]:.0f}" if out_max else "no geometry")
        src_file = find_region_source(hit[1], primary)
        src_max = geojson_max(json.loads(src_file.read_text())) if src_file else None
        if src_max and src_max[0] and out_max:
            factor = out_max[0] / src_max[0]
            check("region scale is a clean pyramid factor",
                  any(abs(factor - c) < 0.01 for c in (0.25, 0.5, 1, 2, 4, 8)), f"x{factor:.2f}")

    if len(cases) > 1:  # synthetic run: prove the panels differ but names resolve
        other = json.loads((projects / "6450.tmap").read_text())
        idx = {layer["name"]: i for i, layer in enumerate(other["layers"])}
        check("batch-2 panel resolves independently",
              idx.get("SST") != names.index("SST") if "SST" in idx else False,
              f"batch-1 SST@{names.index('SST')} vs batch-2 SST@{idx.get('SST')}")

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
