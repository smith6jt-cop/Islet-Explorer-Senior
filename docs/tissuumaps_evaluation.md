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

Those are the real batch-1 and batch-2 layouts from `apps/islet/Channel_names` and
`apps/islet/CLAUDE.md`. **Re-run `scripts/senior/tissuumaps_smoke_test.py` against a real
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

5. **libvips 8.15.1 corrupts reduced pyramid levels when the horizontal predictor is on** —
   and that predictor is its *default* for both `deflate` and `lzw`. Level 0 is always written
   correctly, so the damage is invisible at full zoom and total at every zoom-out. On a sparse
   1%-bright test image whose true mean is 2.00 at every level:

   | `tiffsave` options | level means (libtiff) |
   |---|---|
   | `compression=deflate` (default predictor) | 2.00, **9.65**, **32.69** |
   | `compression=lzw` (default predictor) | 2.00, **9.65**, **32.69** |
   | `compression=deflate, predictor=none` | 2.00, 2.00, 2.07 |
   | `compression=none` | 2.00, 2.00, 2.07 |

   On a real DAPI channel the written levels read 28.0 / 77.0 / 96.0 / 112.5 / 123.7 against
   true block means of 28.0 / 27.5 / 27.0 / 26.0 / 24.9. libtiff and tifffile decode those
   values identically, so the pixels really are wrong in the file — this is not a reader bug,
   and OpenSlide serves the corrupt levels faithfully.

   Visually: background lifts off black into a bright wash as you zoom out, and because
   minority populations lift faster than majority ones the **hue shifts too** — an islet reads
   INS-pink up close and SST-yellow from far out. Composite background measured across levels
   went `(6,5,6) → (168,111,117) → (229,163,164)`, collapsing islet-vs-background contrast to
   nothing. With the fix it stays `(6,5,6) → (7,7,6) → (11,11,11)` and contrast *holds*.

   `predictor="none"` is not a tradeoff: the files also come out **smaller** (DAPI 9.5 → 8.2 MB,
   INS 4.3 → 3.8 MB), because horizontal differencing costs entropy on speckled fluorescence
   rather than saving it. `tissuumaps_smoke_test.py` asserts every reduced level's mean stays
   within 35% of level 0's, which fails loudly without the fix (DAPI 5.2 → 9.7).

## 2a. The white point, and why most markers came out blown

The first robust stretch measured its **white** point over the same population as its **black**
point — everything above `mode + 3σ`. That is right for a black point (background must render
black) but far too low to define "signal": the background's own upper tail, autofluorescence and
haze, reaches well past it.

On a channel where the marker occupied 0.13% of pixels, the population above the black point was
**12.7%** of the image — almost all of it background tail. The marker was therefore the top ~1% of
that set, the 99th percentile landed *below* it, and every marker-positive pixel clipped to flat
white. Fraction of true marker pixels rendered as pure white, scored against ground-truth masks:

| white-point floor | DAPI | INS | GCG | SST |
|---|---|---|---|---|
| `mode + 3σ` (the black point — **was**) | 3.7% | **94.8%** | **90.7%** | 21.2% |
| `mode + 10σ` | 2.1% | 2.8% | 4.6% | 0.6% |
| **`mode + 20σ`** (**now**, `--signal-floor-k`) | 2.1% | 1.4% | 0.9% | 0.2% |
| `mode + 40σ` | 2.1% | 1.3% | 0.9% | 0.2% |

20σ and 40σ agree on every channel: past roughly 20σ the background tail is gone and the estimate
stops moving. The default sits at the point of convergence rather than at a value chosen because
the picture looked nice. This floor is a **display** choice about which pixels inform the white
point — it is not a positivity threshold, and nothing downstream reads it.

**Why the old logs looked fine.** The reported `clip_hi` is a fraction of the *image*, so it read
0.127% on the channel whose every positive pixel had blown out — 0.127% of the image was all the
marker ever was. `clip_signal` now reports the fraction of *unambiguously bright* pixels (fixed at
`mode + 40σ`, deliberately not the floor used to pick the white point, or it would be tautological)
that render as white:

| channel | `clip_hi` (image) | `clip_signal` (marker) — old floor | new floor |
|---|---|---|---|
| INS | 0.130% | **72.7%** | 1.0% |
| GCG | 0.127% | **100.0%** | 1.0% |
| SST | 0.114% | **100.0%** | 1.0% |

**`--report`** prints this table for every channel of a slide from one histogram pass each, writing
no pyramids — so stretch settings can be compared on real data without a re-extraction:

```bash
python scripts/senior/build_tissuumaps_project.py --image /data/islet_ome_tiff/6539.ome.tiff \
    --out-dir /data/tissuumaps --report
```

The smoke test asserts no channel clips more than 5% of its bright population, **and** that the
same check fires at the old floor — a one-sided version would also pass on a fixture too easy to
fail on, which is what the synthetic slides were until they grew an asymmetric background tail and
rare bright markers.

## 2b. Merging mode and layer colours

**Merging mode: use `lighter`.** The Viewer's "Merging mode" dropdown offers 25 options with no
TissUUmaps documentation, because they are not TissUUmaps concepts — `filterUtils.js` passes the
string straight to OpenSeadragon, which sets it as the canvas
[`globalCompositeOperation`](https://developer.mozilla.org/en-US/docs/Web/API/CanvasRenderingContext2D/globalCompositeOperation).
MDN and the W3C Compositing spec are the reference.

QuPath composites fluorescence channels **additively**: each channel's RGB contribution is
summed and clamped at 255. Measured in Chromium against a clamped-sum reference, with two
tinted channels plus a shared dim background:

| mode | result | |
|---|---|---|
| *QuPath additive reference* | `[230, 30, 150]` | |
| `source-over` | `[30, 30, 30]` | TissUUmaps' own default — **top layer only** |
| **`lighter`** | `[230, 30, 150]` | **matches QuPath** |
| `screen` | `[206, 30, 136]` | close, but compresses highlights |
| `multiply` / `darken` | `[0, 0, 0]` | everything vanishes on a dark field |
| `overlay` | `[61, 0, 0]` | |
| `color-dodge` | `[227, 0, 0]` | |
| `difference` | `[170, 30, 90]` | |
| `luminosity` | `[61, 17, 17]` | |

`lighter` is the only exact match, and it is what the generator writes. Note `lighten` (a
per-component max, the other plausible candidate) is **not offered** by the dropdown at all.
Reproduce with the canvas snippet in this section's commit; the modes are browser behaviour, not
ours, so re-measure rather than assume if a browser changes.

**Per-layer sliders: Color, Brightness, Gamma — not Contrast.** TissUUmaps' "Contrast" is not
contrast. OpenSeadragon implements it as `pixel = pixel * adjustment` (`CONTRAST` in
`openseadragon-filtering.js`) — a plain gain with no pivot around mid-grey — on a 0–8 slider. At 2
everything above 128 is already white; at 8, everything above 32. The whole usable range is about
1.0–2.0, which is why small movements look violent, and nothing in it can rescue an over-exposed
layer. `Gamma` is offered instead: `(i/255)^a · 255`, monotonic, clipping nothing, on a 0–2 slider,
and present in the schema enum. Below 1 it lifts dim structure; above 1 it suppresses haze.

**Layer colours are hashed from the marker name.** Only the four default-on markers have fixed
colours (INS red, GCG blue, SST yellow, DAPI dim grey, matching the Avivator palette). Every
other channel gets a fully-saturated hue derived from `md5(marker_name)`.

The previous default was `100,100,100` — a no-op, i.e. **white** — which under additive blending
is as bad as black: it lands on all three components, so switching on a fifth marker both washed
out everything else and looked identical to switching on any other, forcing a manual colour pick
every time.

Hashing the *name* rather than drawing at random or using the channel index is deliberate:

- a **random** colour differs on every rebuild, so a saved figure can't be matched to the live
  viewer, and nothing stops it landing on near-black — the invisible-layer failure fixed in §2;
- an **index-based** colour is worse than random for this cohort, because the panels are ordered
  differently: the same marker would render in different colours on different donors.

Name-hashing gives one stable colour per marker across every donor, panel and rebuild — asserted
by the smoke test, which is the only place two panels exist side by side. The picker's 0–100
integer scale affords roughly 600 fully-saturated hues, so a 59-marker panel sees ~3 collisions
by birthday effect; the builder logs which markers share a colour rather than leaving it to be
discovered, and they only matter when both are displayed together.

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
- **Conversion time.** One-time, offline, idempotent, resumable. Reuse is keyed on an
  `extraction.json` manifest recording the settings *and* an `EXTRACTION_VERSION` bumped
  whenever a code change alters output pixels — not on file mtime alone, which cannot see that
  the code changed and twice let superseded output be served indefinitely. A channel whose
  recorded key no longer matches is re-extracted automatically, with the reason logged.
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
    --regions data/Islets \
    --markers data/app_data/cells/6539_tissue.csv
```

Writes `/data/tissuumaps/6539.tmap`, `channels/6539/<idx>_<marker>.tif`,
`regions/6539.geojson`, `markers/6539.csv`.

`--regions` takes the **directory** the QuPath exports live in (`data/Islets`) and finds the
donor's file by case id, or a single file if you'd rather be explicit. `.gz` is read
transparently, since exports are often stored compressed.

Two rules keep that lookup from attaching the wrong donor's islets, which would land as a
plausible-looking overlay rather than an error. The case id must not be flanked by further
digits — case ids are bare numbers and nothing stops one being a prefix of another, so a plain
`*115*` glob would match donor 1150. And if more than one file still matches, **none** is
attached: the run says so and tells you to pass `--regions <file>`.

**Coordinate conventions** (get these wrong and everything is subtly offset — they are CLI
flags precisely so they can be checked per donor):

| input | space | flag |
|---|---|---|
| viewer OME-TIFF | full-res px, 0.2539 µm/px | `--um-per-px` |
| QuPath GeoJSON | resolution #1 px (`PIXEL_SIZE_UM` = 0.5078) | `--region-um-per-px` |
| per-cell CSV | µm | `coord_factor = 1/--um-per-px`, set automatically |

The region scale factor is **derived, then falsified**, not defaulted. `--region-scale auto`
(the default) computes it as `--region-um-per-px / --um-per-px` — 2.0 for these exports — which
is a physical statement about the two pixel sizes rather than a magic constant. It then checks
the scaled geometry against the image: an outline can never legitimately fall outside the
slide, so a factor that pushes one past the edge is *provably* wrong and the builder backs off
to the largest smaller pyramid factor that fits, logging loudly.

The reverse error — coordinates that are already full-res being doubled — is only detectable
when the tissue covers more than about half the slide. It usually does here, but that is the
known blind spot: the chosen factor and the resulting slide coverage are both printed, and
passing a number to `--region-scale` overrides the inference entirely. The smoke test asserts
that the written outlines land inside the image and that the applied factor is a clean pyramid
factor.

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
# Synthetic: builds both panels, 30 checks, ~2 min. Currently 30/30.
python scripts/senior/tissuumaps_smoke_test.py

# A real donor. --channels keeps the first pass quick; drop it for the full panel.
python scripts/senior/tissuumaps_smoke_test.py \
    --image /data/islet_ome_tiff/6539.ome.tiff \
    --regions data/Islets \
    --markers data/app_data/cells/6539_tissue.csv \
    --channels DAPI,INS,SST,GCG --keep

Rscript scripts/test_viewer_helpers.R                  # 14 unit tests for the R helpers
python scripts/test_build_tissuumaps_regions.py        # 26 unit tests for region handling
```

`test_build_tissuumaps_regions.py` needs no slide, server or pyvips — it covers the directory
lookup, gzip reading, and the scale inference (including the cases the smoke test's single
donor cannot reach: a full-res export, a 4×-downsampled one, and geometry-free input).

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
