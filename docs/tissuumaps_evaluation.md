# TissUUmaps as a replacement for the Avivator viewer

Evaluation + integration notes, July 2026. Question asked: *can the Viewer tab's Avivator
instance be replaced with the more annotation-friendly
[TissUUmaps](https://github.com/smith6jt-cop/TissUUmaps) fork?*

**Answer: yes.** The Viewer tab is already just an `<iframe>` behind a URL builder, so the
app-side change is small. The work is in the image-serving layer. This document records what
was verified, what broke, and how the integration is wired.

---

## 1. What was verified, and on what data

Everything below was exercised on **synthetic slides**, not the real cohort — `data/` and
`/data/islet_ome_tiff/` are not present in a fresh clone (both gitignored, ~1.2 TB of
OME-TIFFs). The synthetic slides were built to reproduce the properties that actually matter
here: **two different panels**, real marker names in the OME-XML, pyramids as SubIFDs (as
`raw2ometiff` writes them), and 16-bit data.

| | synthetic "6539" | synthetic "6450" |
|---|---|---|
| channels | 59 | 58 (no CD38) |
| DAPI / INS / SST / GCG index | 0 / 13 / 35 / 46 | 0 / 14 / **1** / 45 |

Those are the real batch-1 and batch-2 layouts from `app/shiny_app/Channel_names` and
`app/shiny_app/CLAUDE.md`. **Re-run `scripts/senior/tissuumaps_smoke_test.py` against a real
donor slide before trusting any of the size/time numbers below.**

The one thing since confirmed on real data: the generator run against
`/data/islet_ome_tiff/6539.ome.tiff` read 59 channels and resolved `INS@13, GCG@46, SST@35,
DAPI@0` from that slide's own OME-XML — matching the synthetic batch-1 layout exactly. See § 7
for what that run did *not* reach.

Verified end to end (TissUUmaps 3.2.1.14 from the fork, libvips 8.15.1, Chromium headless):

- 59/58 channels detected and extracted to per-channel pyramids; both panels resolve
  DAPI/INS/SST/GCG **by name**, so the default-on channels are correct regardless of panel.
- The generated `.tmap` validates against `tissuumaps-schema` 1.3.
- The server serves every route the project needs: `.tmap` → 200, `.dzi` → 200,
  `_files/<l>/<c>_<r>.jpeg` tiles → 200, markers `.csv` → 200, regions `.geojson` → 200.
- In a real browser: 59 layers named after their markers, exactly DAPI/INS/SST/GCG visible,
  additive (`lighter`) blending, a correct µm scale bar, **6 islet outlines loaded from QuPath
  GeoJSON**, and **2,400 single cells overlaid and grouped by phenotype** with per-phenotype
  counts and colours.

## 2. What broke (and why the app does not use the fork's `/slide` route)

The fork's December 2025 commits add multi-channel TIFF/qptiff support to the `/slide` route.
It works, but it cannot be used as-is for this cohort:

1. **Layer names lose marker identity.** `_extract_multichannel_layers()`
   (`tissuumaps/views.py`) names every layer `<basename>_Channel_<i>` and discards the channel
   names that `_is_multichannel_tiff()` collected. That collection is itself broken: its
   `re.search(r'Name="([^"]+)"', desc)` matches the OME-XML `<Image Name=...>` attribute before
   any `<Channel Name=...>`, so channel 0 came back as `Image0` and the rest as `Channel_1…58`.
   Observed directly:

   ```
   6539: multichannel=True n_channels=59
      channel_names[:6] = ['Image0', 'Channel_1', 'Channel_2', ...]
   layers → [{'name': '6539_Channel_0', ...}, ...]     names carry marker identity? False
   ```

   Positional names are *actively wrong* here, because the cohort mixes panels — SST is index
   35 on batch-1 and index 1 on batch-2. This is the same trap `build_channel_config_b64()`
   avoids on the Avivator side by resolving markers by name.

2. **Extraction runs inside the HTTP request.** `/slide` calls
   `_extract_multichannel_layers()` synchronously, so the first view of a 35 GB slide would
   block a worker for the whole conversion.

Both are avoided by pre-generating the pyramids offline and emitting a `.tmap` that names the
layers — which is what `scripts/senior/build_tissuumaps_project.py` does.

**We run the fork anyway.** The fork's changes are confined to `_fnfilter`, the `/slide` route
and `flask.standalone.addLayer` in `flask.js`, none of which the `.tmap` path touches, so
upstream and the fork serve these projects identically — the fork is not *required*. It is
what we deploy because the viewer should be code we control: any fix we make there (starting
with the two naming bugs above, if we ever want `/slide` to work) ships without waiting on
upstream. `deploy/build_tissuumaps_image.sh` builds the image from
`smith6jt-cop/TissUUmaps`, and `deploy/islet-tissuumaps.container` runs that image.

Two further issues found while building the project generator, both fixed in it:

3. **Regions silently fail to auto-load with >1 layer.** `regionUtils.geoJSON2regions()` opens
   a blocking *"Multiple layers are opened. Please provide the layer index"* modal whenever
   `tmapp.layers.length > 1`, which defeats `autoLoad` — regions stayed at 0 with no console
   error. Stamping `properties.collectionIndex` on every feature skips the prompt entirely
   (`copy_regions()` does this). All channel layers are pixel-aligned, so index 0 is fine.

4. **Greyscale JPEG-in-TIFF is unreliable on libvips 8.15.1.** Single-band pyramid saves fail
   with `VipsJpeg: Unsupported color conversion request` unless the image is built through
   `.scaleimage()` *and* saved with `properties=True`. `deflate` works in every combination, so
   the generator defaults to `--compression deflate` (lossless) and raises an actionable error
   if a `jpeg` run hits the libvips bug.

## 3. Cost of the switch

- **Bit depth and tone mapping.** TissUUmaps layers are 8-bit; Viv renders the 16-bit data with
  live contrast windows. Since marker intensity is not a biological readout for this cohort —
  essentially all of its variation is technical (saturation, background level, and markers whose
  positive fraction spans orders of magnitude) — the generator's job is legibility, not fidelity.
  It therefore defaults to a **background-anchored** stretch (`--stretch robust`): black point at
  the background mode + k·σ (σ from a one-sided MAD below the mode, since signal only contaminates
  the upper side), white point at a quantile of the above-background population with the
  saturation bin excluded, and a channel with nothing above background left dark rather than
  having its noise amplified. A fixed percentile pair cannot do this — it cannot serve a marker
  in 0.05% of pixels and one in 40%. Scored against ground-truth masks (background → signal, as
  rendered 0-255, and how much signal saturates):

  | channel | percentile 0.5/99.5 | robust |
  |---|---|---|
  | dense | bg 3 → sig 160, 1.5% clipped | bg 0 → sig 153, 1.0% clipped |
  | sparse | bg **127** → sig 255, **100% clipped** | bg 0 → sig 177, 3.7% clipped |
  | high background | bg 23 → sig 214, 25.4% clipped | bg 0 → sig 154, 1.3% clipped |
  | saturated | bg 1 → sig 151, 25.2% clipped | bg 0 → sig 121, 8.9% clipped |
  | absent marker | bg **127** — noise rendered as signal | bg 0 — stays black |

  `--stretch percentile` keeps the old behaviour for comparison. Either way this is a display
  proxy — **not** a substitute for the quantitative pipeline.
- **Disk.** A second copy of every channel. On synthetic data the ratio is meaningless (smooth
  gradients compress absurdly well); the smoke test prints the real ratio when pointed at a
  donor slide. Budget for it before converting all 22.
- **Conversion time.** One-time, offline, idempotent (skip-if-exists), resumable.
- **A second service.** Trades the nginx sidecar (`islet-images.container`) for the TissUUmaps
  container. Net zero services — and the HTTP-Range requirement disappears, so the Viewer works
  under `shiny::runApp` locally, which it currently does not.
- **Security.** The Flask app exposes a file browser over its slide directory. The deployment
  sets `READ_ONLY = True` and mounts the volume `:ro`.

## 4. What you gain

- Region drawing (free-hand, polygon, brush, rect, ellipse), rename/reclassify, boolean ops,
  and **GeoJSON import/export** — so QuPath islet boundaries render *on the image*, and edits
  come back out as GeoJSON.
- Per-cell marker overlay straight from the existing CSVs (`X_centroid`, `Y_centroid`,
  `phenotype`), millions of points, GPU-rendered, grouped and coloured by phenotype.
- Region statistics — area, perimeter, and marker counts inside any drawn ROI, exportable as
  CSV. That is per-islet composition on arbitrary user-drawn regions, which the drill-down can
  only do on fixed QuPath annotations today.
- A route to zoom-to-islet: the `.tmap` schema has a `boundingBox` key, so a future
  `forced_image` equivalent could deep-link to a clicked islet (Avivator only ever selected the
  file).

## 5. How it is wired

**Install the fork.** Everything below assumes the TissUUmaps you run is
[`smith6jt-cop/TissUUmaps`](https://github.com/smith6jt-cop/TissUUmaps), not the PyPI release
or the published `cavenel/tissuumaps` image.

```bash
# Local / dev — into the same env that runs the project builder.
# Needs libvips + openslide present (pyvips and openslide-python are bindings, not
# self-contained): apt install libvips libopenslide0, or
# conda install -c conda-forge libvips pyvips openslide-python.
pip install "git+https://github.com/smith6jt-cop/TissUUmaps.git@master"
tissuumaps_server -l 127.0.0.1 -p 5100 /data/tissuumaps

# Production — build the container image from the fork (~10-20 min cold; the fork's
# Dockerfile compiles openslide 3.4.1 from source). Tags :<short-sha> and :latest.
bash deploy/build_tissuumaps_image.sh
```

`tissuumaps_smoke_test.py` prints the version and install path of whatever TissUUmaps it
exercises, so a run against a pip-installed upstream is visible in the log rather than assumed.

**Build a project per donor** (idempotent; re-run after any slide or GeoJSON change):

```bash
python scripts/senior/build_tissuumaps_project.py \
    --image   /data/islet_ome_tiff/6539.ome.tiff \
    --out-dir /data/tissuumaps \
    --regions data/app_data/json/6539.geojson \
    --markers data/app_data/cells/6539_tissue.csv
```

Writes `/data/tissuumaps/6539.tmap`, `channels/6539/<idx>_<marker>.tif`,
`regions/6539.geojson`, `markers/6539.csv`.

**Coordinate conventions** (get these wrong and everything is subtly offset — they are CLI
flags precisely so they can be checked per donor):

| input | space | flag |
|---|---|---|
| viewer OME-TIFF | full-res px, 0.2539 µm/px | `--um-per-px` |
| QuPath GeoJSON | resolution #1 px (`PIXEL_SIZE_UM` = 0.5078) | `--region-scale 2.0` |
| per-cell CSV | µm | `coord_factor = 1/--um-per-px`, set automatically |

**Point the app at it** — `ISLET_VIEWER_BACKEND=tissuumaps` plus `TISSUUMAPS_URL` and
`TISSUUMAPS_PROJECT_DIR` (see `deploy/env.sh.example`). Anything other than `tissuumaps` keeps
Avivator, so this is a per-deployment flip, not a migration.

**Serve it** — `deploy/islet-tissuumaps.container` (runs the fork-built
`localhost/islet-tissuumaps:latest`) + `deploy/tissuumaps.cfg`, and a `/tissuumaps/` proxy
route alongside the existing app route. Because the image is built locally rather than pulled,
it must exist on the host before the unit starts, and picking up new fork commits means
re-running the build script and restarting the service — there is no auto-update.

## 6. Re-running the verification

```bash
# Synthetic: builds both panels, 19 checks, ~2 min. Currently 19/19.
python scripts/senior/tissuumaps_smoke_test.py

# A real donor. --channels keeps the first pass quick; drop it for the full panel.
python scripts/senior/tissuumaps_smoke_test.py \
    --image /data/islet_ome_tiff/6539.ome.tiff \
    --regions data/app_data/json/6539.geojson \
    --markers data/app_data/cells/6539_tissue.csv \
    --channels DAPI,INS,SST,GCG --keep

Rscript scripts/test_viewer_helpers.R      # 14 unit tests for the R helpers
```

The smoke test builds (or reuses) a slide, runs the generator, boots `tissuumaps_server`,
asserts every route, and drives a headless Chromium to confirm layers/regions/markers actually
render. `--keep` leaves the server up so you can look at it. `SMOKE_SHOT=/path/shot.png` saves
a screenshot; `CHROMIUM_PATH` overrides the browser binary.

One environment gotcha seen in a Debian/Ubuntu container (not conda): importing `pyvips` pulls
in the system `libhdf5` 1.10, which then breaks `h5py`'s bundled 2.0 with
`ValueError: Not a datatype`. Workaround if you hit it:
`LD_PRELOAD=$(python -c "import h5py,glob,os;print(glob.glob(os.path.join(os.path.dirname(h5py.__file__),'..','h5py.libs','libhdf5-*.so*'))[0])")`.

## 7. Not done

- **Channel-name resolution is confirmed on a real slide; nothing else is.** A first run of
  the generator against `/data/islet_ome_tiff/6539.ome.tiff` read 59 channels and resolved
  `INS@13, GCG@46, SST@35, DAPI@0` — the batch-1 layout, straight from that slide's OME-XML.
  The run then stopped on a missing `pyvips`, so **no real channel pyramid has been extracted
  yet**: every size, ratio and timing figure in this document is still synthetic, and the 8-bit
  rendering has not been eyeballed on tissue.
- The container image was not exercised (no Docker daemon in the dev container); verification
  used `tissuumaps_server` from the fork's source, which is the same Flask app the image runs
  under gunicorn.
- Avivator is still the default and nothing has been deleted. Removing
  `www/avivator/`, `scripts/install_avivator.sh`, `scripts/make_offsets_json.py`, the
  `.offsets.json` sidecars and the nginx sidecar is a follow-up, once the 8-bit rendering has
  been eyeballed on real tissue.
- Zoom-to-islet (`boundingBox` deep-linking) is not implemented.
