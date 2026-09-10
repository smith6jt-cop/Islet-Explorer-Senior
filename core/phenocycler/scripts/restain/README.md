# Post-Xenium IF restain lane

## Current registration-only cohort workflow

This lane registers each post-Xenium IF qptiff to the matching pancreas Xenium morphology grid and
publishes an already-warped OME-TIFF for Xenium Explorer's **Yes, skip alignment** route. When a
same-section QuPath GeoJSON exists, it also transforms the filtered islet outlines to an
Explorer-importable annotation sidecar on that exact grid. The post-Xenium IF images the same
physical section Xenium measured, so Xenium's nucleus polygons and the IF DAPI peaks describe the
same nuclei. No new cell segmentation, cell typing, positivity calling, or biomarker
interpretation is part of this workflow.

Inputs and working outputs live outside the content-addressed Phenocycler pipeline. Generated
state is under `data/development/restain/`; this lane does not write a `StageManifest`, touch
`data/runs/`, or contribute scripts under `scripts/**` to a `run_id`.

### Input and cohort contract

`data/integration/restain_registration_inputs.csv` is the source inventory for registration. It
maps every transferred qptiff to an exact donor, pancreas Xenium serial, ROI, and provenance
source. Before work starts, `registration_inputs.py` checks that:

- the inventory equals the local qptiff set, with one image, donor, and serial per row;
- the serial resolves to the same donor in `data/integration/xenium_paths.csv`;
- the qptiff has one literal DAPI page, `uint8` data, a floor-halved pyramid, named raw channels,
  and one exact XML `PixelSizeMicrons` value; and
- the matching Xenium bundle has `cells.parquet`, `nucleus_boundaries.parquet`, and the DAPI
  morphology file.

The inventory deliberately contains acquisition records that are not executable cohort members.
Donors `6457` and `6579` remain centrally excluded by `phenocycler/cohort.py`: an explicit request
for either donor fails, and an omitted `--donors` filters them out. Local stale VALIS directories
cannot bypass that guard.

QuPath GeoJSON is optional for image registration. `nuclei.py` uses a drawn `Tissue` polygon when
that class exists, records `whole_slide_no_geojson` when no file exists, and records
`whole_slide_no_tissue_annotation` for an islet-only file. `run_align.py --windows auto` likewise
uses tissue windows only when the file actually contains `Tissue`; otherwise it retains the
annotation-independent occupied-window gate. `--whole-slide` can force nucleus detection to ignore
annotations. Thus an islet-only export is downstream geometry, not a new registration input, and
cannot perturb the already-validated DAPI transform. DAPI export, VALIS, nucleus refinement,
occupied-window gating, and registered OME export do not require antibody roles. The five sections in
`restain_sections.csv` remain the only sections with maintainer-declared biological channel roles;
the registration inventory does not extend those declarations.

`registration_inputs.py` pairs a GeoJSON only by the same Xenium serial as the qptiff. Annotation
vertices are interpreted as qptiff level-0 pixel edges and scaled with that qptiff's exact XML
`PixelSizeMicrons`; there is no hard-coded annotation scale in the export path.

### DAPI seed and nucleus refinement

Run from `core/phenocycler/`. With `--donors` omitted, each command selects every centrally
eligible inventory row. For a resumable cohort run, pass the same donor subset to each command;
the two-donor subset below is only command syntax, not a claim that either result passes.

```bash
cd core/phenocycler
../../.venv/bin/python scripts/restain/export_dapi.py --donors 115 6374
~/.venvs/valis/bin/python scripts/restain/valis_register.py --donors 115 6374
../../.venv/bin/python scripts/restain/nuclei.py --donors 115 6374
../../.venv/bin/python scripts/restain/run_refine.py --donors 115 6374
```

`export_dapi.py` selects pyramid levels from each image's exact calibration and records the two
different pixel-origin conventions in `export.json`. VALIS supplies a DAPI-to-DAPI seed in its
isolated environment; failure of its optional micro pass is recorded while retaining the rigid
seed for downstream measurement. `nuclei.py` detects IF DAPI peaks and reads Xenium nucleus
centroids from `nucleus_boundaries.parquet`. `run_refine.py` reports every candidate and spatial
hold-out, then writes the chosen `refined` transform. A small in-sample nearest-neighbour residual
is not a release decision.

### All-window gate

Gate every candidate on the annotation-free, origin-anchored 500 µm occupied-window grid. Specify
`occupied` explicitly so annotated and unannotated donors use the same window contract, and use
all windows before deciding whether an intensity correction is warranted:

```bash
../../.venv/bin/python scripts/restain/run_align.py \
  --donors 115 6374 --transform refined --windows occupied --window-fold all
```

The gate rasterises transformed Xenium nucleus polygons and measures their DAPI displacement by
normalised cross-correlation in every occupied window. It writes a donor-local
`alignment_refined.json` even for `FAIL` or `INCONCLUSIVE`, including the reason, measured and
unmeasured windows, bounded residual, contrast checks, transform hashes, and thresholds. Exit
codes are `0` for `PASS`, `3` for `INCONCLUSIVE`, `4` for `FAIL`, and `1` for an I/O/setup error.
A donor that does not pass has no releasable registered image.

### Optional one-fold intensity refinement

Use this only when the all-window `refined` gate exposes a coherent DAPI offset field. The fitter
reads one deterministic checkerboard fold, fits one robust affine residual represented on a coarse
field, and never reads the other fold. The example trains on `even`, so the release gate must use
the held-out `odd` fold:

```bash
../../.venv/bin/python scripts/restain/run_intensity_refine.py \
  --donors 115 --source-stem refined --output-stem intensity_refined --training-fold even
../../.venv/bin/python scripts/restain/run_align.py \
  --donors 115 --transform intensity_refined --windows occupied --window-fold odd
```

The correction is optional, refuses insufficient or implausible training evidence, writes a new
`intensity_refined` transform rather than overwriting `refined`, and records its training/withheld
fold contract. Reversing the training fold requires gating `even` instead. The registered exporter
rejects an intensity-refined transform unless `alignment_intensity_refined.json` is `PASS` on the
declared held-out fold; a full-window or training-fold PASS is not interchangeable with that test.

### Exact-stem registered OME export

Export donors in separate calls grouped by the transform that passed. `--transform-stem` reads
exactly `<stem>.json`, its optional `<stem>_displacement.npz`, and
`alignment_<stem>.json`; it never falls back to another transform or reconstructs a transform from
the gate. Donor/serial/input identity, direction, calibration, affine, displacement presence, and
content hashes are checked before the first multi-gigabyte write.

```bash
# Donors whose all-window nucleus-refined gate passed
../../.venv/bin/python scripts/restain/run_registered_export.py \
  --donors 6374 --transform-stem refined

# Donors whose optional correction passed on its withheld fold
../../.venv/bin/python scripts/restain/run_registered_export.py \
  --donors 115 --transform-stem intensity_refined
```

Each successful publication is
`data/development/restain/explorer_bundle/donor_<id>_<serial>/`, containing:

- `images/<serial>_IF_on_xenium_grid.ome.tif`, with every raw qptiff channel in native page
  order, resampled through the full selected transform to the exact Xenium morphology dimensions,
  `0.2125 µm/px` scale, and pyramid level shapes;
- `manifest.json`, which binds the image to the input, PASS gate, exact transform, channel order,
  and import route;
- `qc.json`, which verifies OME metadata, BigTIFF, shape, dtype, scale, tiling, compression, and a
  decoded tile from every channel at every level; and
- a donor-local `README.md` with the Xenium Explorer import instruction.

When the inventory finds a same-serial GeoJSON, the bundle additionally contains
`annotations/<serial>_islets_on_xenium_grid.geojson`. Source QuPath pixels are converted to IF
microns, every exterior vertex is mapped through `Transform.inverse_apply()` (the stored transform
maps Xenium microns to IF microns), and the result is divided by `0.2125 µm/px` to land on the
Xenium morphology level-0 pixel grid. The full inverse includes a displacement field when the
selected transform has one.

Only `Islet` polygons at least 50 µm in equivalent diameter are emitted; sub-cell classifier
fragments, other classes, and unsupported source geometry are counted in the manifest. Before
publication, the file is reopened and checked against Xenium Explorer 4.1's third-party GeoJSON
contract: a nonempty `FeatureCollection` of fewer than 10,000 closed, simple, hole-free `Polygon`
features, fewer than 100,000 vertices each, finite coordinates inside the morphology grid,
`objectType: annotation`, and an `Islet` classification with RGB color. A malformed sidecar aborts
the transaction before the image render begins.

In Xenium Explorer 4.1 or newer, open the matching Xenium bundle, import the OME-TIFF with **Yes,
skip alignment**, then use **Annotations → Import** for the GeoJSON. Associate it with the **Xenium
In Situ morphology image** and choose **pixels** if asked for units. Never apply an image-alignment
CSV to either already-transformed product.

Publication is transactional through a sibling `.partial` directory. Existing bundles are
refused by default; `--skip-existing` revalidates and keeps one, while `--overwrite` replaces a
non-protected bundle only after its replacement passes QC. Donors `6436` and `6442` are protected
completed Explorer exports and are never rewritten by this driver, even with `--overwrite`.

To add newly arrived annotations to an existing bundle without regenerating its multi-gigabyte
OME-TIFF, use the annotation-only publisher:

```bash
../../.venv/bin/python scripts/restain/run_islet_annotation_export.py --donors 6436 6442
```

It verifies bundle, qptiff, morphology, and serialized-transform identity; builds and validates a
sibling partial sidecar; then publishes only the annotation metadata/files. The two historical
protected bundles are eligible only through this narrow path, and only when their manifest's exact
transform still equals the donor-local serialized transform. Existing annotation sidecars are
refused unless `--overwrite-annotations` is given; images and cell groups remain untouched.

The qptiffs expose technical names such as `DAPI`, `Opal570`, `Opal620`/`Opal690`, `Opal520`, and
`Sample AF`, but no reliable Biomarker metadata. The registered OME preserves those raw names and
does **not** infer INS, GCG, or another biological identity. Its manifest explicitly leaves every
role unassigned and contains no positivity calls or `cell_groups` directory.

### Per-donor audit and failure records

The resumable authority is the donor directory
`data/development/restain/valis/donor_id=<donor>/`, not a batch summary that a later subset can
replace.

| file | durable evidence |
|---|---|
| `export.json` | exact qptiff/Xenium inputs, DAPI pages, levels, calibration, pixel origins, and image extents |
| `valis_result.json` | seed provenance, VALIS error summary, optional-micro status, and warped-centroid displacement |
| `nuclei.json` | IF/Xenium nucleus counts, detector settings, annotation/whole-slide policy, and density ratio |
| `refinement.json` | every point-set candidate, selected winner, detection ceiling, placement null, and spatial hold-out |
| `refined.json` plus optional `refined_displacement.npz` | exact nucleus-refined transform serialization |
| `alignment_<stem>.json` | per-window release gate, exact transform hashes, `PASS`/`INCONCLUSIVE`/`FAIL`, and `verdict_why` |
| `intensity_refinement_<stem>.json` | optional correction fit, source hashes, training fold, withheld fold, and required gate command |

A fatal missing-input or malformed-file error can stop before a stage has evidence to write; it is
not converted into a synthetic PASS/FAIL record. Conversely, a measured registration failure stays
on disk as the donor's gate report and the exporter refuses it. Do not force, relabel, or omit that
failure to complete a nominal cohort count.

## PhenoCycler annotations on the Xenium grid (serial section)

The sections above are all same-slide: the post-Xenium IF images the section Xenium measured, so
the correspondence is a measurement. The 55-plex PhenoCycler section is **serial** to Xenium, and
its drawn anatomy reaches Explorer by a different route — the S4 registration in
`phenocycler/integration/register.py`, not the nucleus-matched transform above. The output shares
the Explorer 4.1 annotation contract and nothing else about what it may claim.

```bash
cd core/phenocycler
../../.venv/bin/python scripts/restain/run_pheno_annotation_export.py --donors 6436 6442
```

The driver is annotation-only and additive: it stages, validates every file, then moves them into
the existing bundle and rolls the whole set back on any failure. It never touches the registered
OME-TIFF, alignment files, keypoints, cell groups, or the same-slide IF islet sidecar. Published
files carry `pheno` in the name — `<serial>_pheno_<Class>_on_xenium_grid.geojson` — because a
reviewer must be able to tell a serial-section outline from a same-section one in Explorer's layer
list without opening the QC file. `pheno_annotations_qc.json` records every count and the
registration evidence behind it.

### Three things this lane handles that the same-slide lane does not

**The Xenium section is a portion of the PhenoCycler section.** On the two registered donors the
Xenium extent is 39.4% (6436) and 54.2% (6442) of the PhenoCycler extent, and lies entirely inside
it. So annotations outside the Xenium field have nowhere to land. They are dropped and counted,
never clipped — clipping would invent an anatomical boundary that is not there and would misstate
the annotation's area.

**The field is the tissue outline, not the image rectangle.** Until 2026-09-10 the only spatial
test was the morphology image's bounding rectangle. That rectangle is the capture area, and the
tissue is **72.3%** of it on *both* donors — so more than a quarter of the area an annotation
could pass the test in is blank background Xenium never imaged, and annotations were reaching
Explorer sitting on it. `PA.xenium_tissue_outline` builds the region from the Xenium **cell
centroids** through `rasterize.tissue_mask` — the same function `register.evaluate` scores tissue
Dice with, so the region annotations are selected against and the region the registration is
judged on cannot drift apart — polygonises it with `cv2.findContours`, and the exporter keeps only
polygons lying **entirely** inside it. Note `tissue_mask` fills holes, so an interlobular gap
inside the section is not treated as a place Xenium failed to look. The outline, its provenance
and the per-class counts it removed are all recorded in `pheno_annotations_qc.json`, and the
outline itself is published as its own Explorer layer so a reviewer can see the region every drop
decision was made against.

The extent fact also reframes `qc.py`'s verdict, though not in the way this file used to claim.
Whole-section tissue Dice cannot exceed **0.565** and **0.703** for these donors on extent alone,
so a 0.80 bar was unreachable by construction. But the deeper problem is not the ceiling: the
metric moves *against* correctness — see "Is a donor merely mis-oriented?" — so raising or
lowering the bar could not have fixed it. Tissue Dice is now reported and never gated.

**The non-rigid field is not invertible here.** The B-spline stage is off for this lane as of
2026-09-10 (`reg_nonrigid = false`). Measured on the real annotation vertices under the corrected
transforms, **15.0%** (6436) and **6.1%** (6442) fail the round trip the exporter requires, which
costs **90** and **48** publishable islets; affine-only round-trips at 4e-12 µm on the same
vertices. The earlier fields were worse than non-convergent — they were *folded*: 6436's sat at
its 200 µm cap on 14.7% of pixels with a minimum Jacobian determinant of −1.657 over 0.85% of the
field, so `x -> x + D(x)` was not injective there and the fixed-point inverse had no solution to
find. The per-pixel radial clip in `refine_bspline` is the mechanism, not the safeguard its
docstring claimed. `displacement_report` now records the Jacobian and the round trip, and
`displacement_is_invertible` gates on them rather than on `max_displacement_um`, which cannot see
folding at all.

**One QuPath class exceeds what Explorer will import.** Explorer takes fewer than 10,000 polygons
per file. After dropping what falls outside the Xenium tissue, `Blood_vessel` still carries 76,159
(6436) and 88,578 (6442) parts, and no defensible diameter floor brings that under the cap.
Rather than tune a filter to a display limit, a class over the cap is **split into ordered parts,
largest object first**, and every polygon is published. `--min-diameter-um` remains available,
defaults to 0, and is recorded when used.

### Scale is the thing that has to be right

PhenoCycler and Xenium both report **calibrated microns**. So the transform between two serial
sections of one block is physically rigid, up to whatever real stretch sectioning and mounting
introduced: scale about 1 on *both* axes, negligible shear. Anything else is a fitting artefact,
and it is the one error that hides in a whole-section overlay while corrupting every micron
measurement taken through the transform.

`Transform.scale` reports a single isotropic number — the square root of |det| — which is exactly
where a single-axis stretch hides. Decomposed per axis:

| | 6436 | 6442 |
|---|---:|---:|
| scale x / scale y | **1.1057 / 0.9989** | 1.0020 / 0.9898 |
| anisotropy (sx/sy) | **1.107** | 1.012 |
| isotropic `Transform.scale` | 1.0509 | 0.9958 |
| shear | 0.55° | 0.10° |
| drift over a 10 mm section | **509 µm** | 42 µm |

6436's isotropic scale of 1.0509 looks unremarkable; it is hiding an **11% x-only stretch**.

The source is `refine_ecc(motion="affine")`. ECC's affine mode has six free parameters, and driven
by a partial-overlap silhouette it will buy correlation with a stretch. The preceding
`motion="euclidean"` stage is rotation and translation only, so its scale is exactly 1 by
construction. `scale_is_physical` therefore gates the affine refinement: it is taken only when it
both scores better *and* stays inside `reg_max_anisotropy` (1.03) and `reg_max_scale_deviation`
(0.05); otherwise the rigid fit is kept and the refusal is recorded in the donor's `notes`. Every
donor's `qc.json` now carries a `scale_report` block with both axis scales, so this is auditable
per donor rather than something to rediscover.

This matters more than the last few degrees of orientation. A registration that is a little off in
rotation misplaces an annotation by a bounded amount everywhere; one that is off in scale misplaces
it by an amount that grows across the section, and no downstream consumer can detect that from the
output alone.

### Is a donor merely mis-oriented? Ask the islets — but ask about translation too

**This section is a correction.** It previously concluded that donor 6436's two sections did not
share an islet field, and that conclusion was wrong. 6436's transform was out by 12.5 deg and
2.6 mm; corrected, its islet correspondence goes from 12 of 136 Xenium islets within 100 um —
against a chance expectation of 5 — to **69 of 136**, which is the same fraction 6442 achieves.
The prose that stood here is preserved in Git history; what follows is why it was wrong, because
the reason generalises.

Rotation is the first thing to suspect when a serial-section registration looks off, and
`register.py` does search it exhaustively (0-360 deg at 5 deg, plus the mirror). The diagnostic
built to check it, `pheno_rotation_figures.py`, sweeps a residual rotation on top of the stored
transform and plots islet Dice and islet RMSE — both from islet centroids, which the image track
never sees — beside the footprint-coverage curve. For 6442 those curves spike sharply at 0 deg.
For 6436 they were flat at every angle, and that was read as "no rotation recovers a
correspondence that is not there".

**A rotation-only sweep cannot support that reading.** It rotates about the warped moving
centroid, which fixes that centroid by construction, so it traces a one-parameter curve through a
three-parameter space and no point on it carries a translation. Both metrics on it also saturate:
islet Dice dilates each side by 75 um and so reads zero past ~100 um of residual, and islet RMSE
is capped at 300 um on mutual pairs, so under a wrong transform it collects chance pairs and
converges to the null rather than to a bad value. 176.9 um on 52 pairs *was* the null. A joint
(rotation, translation) error of the size actually present produces exactly the flat curves that
were observed — so the curves were consistent with the conclusion drawn, and equally consistent
with the transform simply being in the wrong place.

`pheno_residual_search.py` asks the question the sweep could not, by making translation free at
every angle:

```bash
cd core/phenocycler
../../.venv/bin/python scripts/restain/pheno_residual_search.py --donors 6436 6442
```

For each residual angle over the full 360 deg x mirror it bins every moving-to-fixed offset
vector and takes the peak, so the best translation for that angle is read off rather than
optimised for. The statistic is `peak / p95` of the donor's own landscape — a within-donor null,
carrying no expectation about how many islets *should* match.

Measured against the **stored transform of the day**:

| | 6442 | 6436 |
|---|---|---|
| peak location | **0 deg, (-25, +25) um** | −14 deg, **(−675, +425) um** |
| peak contrast | 6.96x the landscape p95 | 5.68x |
| mirror | none | none |
| Xenium islets within 100 um | 66 of 143 | **12 of 136** (chance ~5) |

6442 recovering zero residual to within one 50 um bin is what makes 6436's 800 um peak
trustworthy: the same statistic, on the same code, says "already right" for one donor and "800 um
that way" for the other. Ducts (1 -> 4 of 28) and vessels (5 -> 15 of 110) move with it, and
neither is hormone-seeded, so neither is confounded by 6436 being T1D.

**What actually went wrong.** The coarse stage ranked 6436's correct hypothesis **4th of 144**,
and `register_image` refined `candidates[:n_coarse_candidates]` with `n_coarse_candidates = 3`.
It was cut by one position. The ranking statistic was not the problem — it put the truth 1st for
6442 and 4th for 6436 — but its margins are small when one section is a sub-field of the other,
because most placements overlap a similar amount of tissue whichever way the moving image is
turned. So the cut is now a plateau rather than a rank (`reg_coarse_keep_frac`, floored at
`reg_coarse_min_candidates`), and ECC is masked to the Xenium footprint, which both widens the
margin and widens the basin it can converge from.

The earlier overlap-restricted NCC experiment, recorded here as having "made every independent
metric worse", drove 6436 to −30.4 deg — which *crosses* the true −19.6 deg. It was judged by
islet Dice and islet RMSE with no translation re-optimisation, i.e. by metrics structurally
incapable of registering the improvement. The reason it was rejected was not the reason it was
wrong. It is still not the fix — the masked coarse statistic fails the 6442 positive control —
but the record should say why.

**The lesson that generalises: no silhouette metric can grade this.** At 6436's correct transform,
tissue Dice *falls* from 0.549 to 0.530 and footprint containment falls from 0.910 to 0.890. Two
serial-section outlines can be brought into agreement by many placements, and the one that agrees
best is not the one that matches anatomy. `score()` weighted tissue Dice at 0.5 and `qc.py` gated
on it at 0.80; both therefore preferred the wrong answer, and both have been replaced — see
"What the gate is now".

### What the gate is now

Two hard preconditions, then one primary bar and one corroborating bar:

* **Scale physicality** (`scale_is_physical`) — unchanged, and it works: it caught 6436's 11%
  x-only stretch, which was itself a symptom of the wrong basin. With the basin corrected the
  affine is accepted at anisotropy 1.005.
* **Displacement invertibility** (`displacement_report`, `displacement_is_invertible`) —
  `max_displacement_um` cannot see folding, and the per-pixel radial clip in `refine_bspline` is
  what causes it. 6436's first field sat at its 200 um cap on 14.7% of pixels, folded on 0.85%
  (min Jacobian det −1.657), and missed its round trip by up to 400 um. Gate on the Jacobian and
  the round trip, both pure geometry.
* **Primary: `footprint_ncc_z`** — correlation between the two DAPI images inside the warped
  Xenium footprint, scored against 120 other placements of that same footprint on that same
  section. Within-donor, so staining, density and disease state are held fixed and only placement
  varies. Never gate the absolute value: 6436's correct placement reads 0.303 and 6442's reads
  lower, so a bar set from either would reject the other. Measured: z=6.9 (6436) and z=6.7
  (6442), against z=2.3 for 6436's old transform.
* **Corroborating: `islet_peak_ratio` / `islet_peak_offset_um`** — does the joint landscape peak
  *here*. Both donors now peak at 35 um and ~6x p95. Corroborating only, never decisive: both
  islet callers are hormone-seeded and weight a T1D donor's islets differently across the two
  modalities.

`reg_nonrigid` is off for this lane as of 2026-09-10. Not preference — measurement. On the real
annotation vertices, 15.0% (6436) and 6.1% (6442) failed the round trip the exporter requires,
costing 90 and 48 publishable islets; affine-only round-trips at 4e-12 um on the same vertices.

### Where both donors stand now

| | 6436 before | **6436 after** | 6442 before | **6442 after** |
|---|---:|---:|---:|---:|
| rotation | −7.13 deg | **−19.60 deg** | +5.11 deg | **+5.09 deg** |
| scale (sx / sy) | 1.0000 / 1.0000 | **0.9931 / 0.9880** | 1.0020 / 0.9898 | **1.0030 / 0.9902** |
| islet RMSE | 176.9 um | **89.1 um** | 107.9 um | **92.1 um** |
| Xenium islets within 100 um | 12 of 136 | **69 of 136** | 66 of 143 | **66 of 143** |
| ducts within 150 um | 1 of 28 | **4 of 28** | 11 of 27 | **11 of 27** |
| vessels within 150 um | 5 of 110 | **15 of 110** | — | — |
| footprint NCC z | 2.3 | **6.9** | — | **6.7** |
| islet landscape peak | 800 um away | **35 um** | 35 um | **35 um** |
| tissue Dice | 0.549 | 0.530 | 0.767 | 0.768 |
| verdict | FAIL | **PASS** | WARN | **PASS** |

6442 barely moves, which is the regression control: the same code change that relocated 6436 by
2.6 mm left the donor that was already right where it was. Its islet RMSE improving 107.9 -> 92.1
um comes from the masked ECC and the two coarse-stage corrections.

Note that islet Dice is *not* comparable across this change, because the dilation element was
fixed in the same commit — scipy's default is the 4-connected cross, which iterates to an L1
diamond, not the 75 um disc the docstring claimed. Holding the metric fixed, 6436's islet Dice
goes 0.024 -> 0.178 on the transform alone; the element fix roughly doubles it for any transform.

### What was published

Six classes per donor. `Islets` carries two — the core and its 20 µm expansion, matching `Islet_N`
/ `Islet_N_exp20um` upstream. Counts are polygons landing wholly inside the Xenium **tissue
outline**, through the corrected transforms:

| class | 6436 | 6442 |
|---|---:|---:|
| `Islet` | 208 | 293 |
| `IsletExpanded` | 212 | 295 |
| `Tissue` | 411 | 191 |
| `Blood_vessel` | 76,159 (8 files) | 88,578 (9 files) |
| `Lymphatic` | 6,561 | 25,175 (3 files) |
| `Nerve` | 1,475 | 4,113 |
| **total** | **85,026** in 13 files | **118,645** in 16 files |

Two changes moved these counts at once, and the QC file separates them per class. The corrected
transform relocated 6436's field by 2.6 mm, and the tissue-outline test removed what the image
rectangle had let through: 24,867 further parts on 6436 and 24,434 on 6442, all of them
annotations that had been landing on blank background. Republishing retires this lane's own
previous files for the classes being written, so a class that splits into fewer parts than last
time cannot leave a stale part behind carrying a superseded transform.

**What this is not.** A serial-section registration locates anatomy approximately. These outlines
are not a cell-level correspondence, they do not transfer any PhenoCycler identity or measurement
onto a Xenium cell, and they must not be read as one. The islet RMSE behind them is 89 µm (6436)
and 92 µm (6442).

Perfect alignment is not the bar here and is not achievable between serial sections. Accurate
*scale* is a bar, because a scale error is the one kind that grows across the section and that a
downstream consumer cannot detect from the output: 6436 is at 0.9931/0.9880 and 6442 at
1.0030/0.9902, both inside the anisotropy bound. Being in the right *place* is the other bar, and
it is the one this lane got wrong for a month while every scale check passed.

## Historical verified baseline: donors 6436 and 6442

The sections below preserve the biological investigation and the two registrations already
reviewed in Xenium Explorer. Their numbers establish the historical two-donor baseline; they are
not a status table for the remaining cohort and do not imply that every acquisition can pass the
same registration gate.

Per-cell: the distance from a Xenium nucleus to the IF nucleus it lands on, over every nucleus in
the section.

| donor | tissue outlines (ECC) | VALIS, as it stood | **nucleus-matched** | p90 | within 2 µm | null |
|---|---:|---:|---:|---:|---:|---:|
| 6436 | 42.13 µm¹ | 2.84 µm | **0.546 µm** | 1.571 | 53.1 % | 8.1 % |
| 6442 | 31.46 µm | 1.36 µm | **0.630 µm** | 1.653 | 57.1 % | 8.0 % |

¹ window residual; only 25 of 243 windows were even measurable at that alignment.

`within 2 µm` is the share of **all** Xenium nuclei with an IF nucleus that close. Its ceiling is
the IF detector's own recall — 56.7 % and 61.0 % — so 53.1 % and 57.1 % are 94 % of what was
available to match, and the placement null is 8 %.

The gate (`run_align.py`, below) returns **PASS on both donors**, at a median window residual of
**0.11 µm (6436)** and **0.13 µm (6442)** with every window measurable, against a bar of 1.80 and
1.88 µm.

## How, and what was wrong before

**The registration is `register.icp` on the two nucleus fields**, seeded by VALIS. Xenium's side is
`nucleus_boundaries.parquet` — already a segmentation of Xenium's own DAPI — and the IF side is
detected by `nuclei.py`: a band-pass at the ~6 µm nucleus scale, local maxima no closer than 4 µm,
and a parabolic sub-pixel fit. 395,048 and 584,723 IF nuclei inside the drawn tissue. No cell
objects and no cell types are created; these points exist only to be matched.

The seed only has to land inside half the spacing between nuclei (~4.5 µm) for the mutual
nearest-neighbour matching to be unambiguous. Three defects had been keeping VALIS outside that.
The third is pinned by tests; the first two are inside VALIS's own call and cannot be tested from
the pipeline venv, so what stands for them is VALIS's own reported error, recorded below and in
`valis_result.json`:

**1. VALIS's rigid transform had never been refined above a thumbnail.** `Valis.register` fits `M`
on the processed images — 512 px on a 20,880 px slide, about 40 µm/px. The only thing in
valis-wsi 1.2.0 that re-fits `M` is `MicroRigidRegistrar`, which `Valis.__init__` defaults to
`None` (registration.py:1790) and which `register()` runs between the rigid and non-rigid stages
(registration.py:4743-4745). Enabling it: `N low rez matches = 1372, N high rez matches = 31596,
low rez D = 1.190, high rez D = 0.231`, and VALIS's own error fell

    6436   rigid_D 15.28 -> 3.892 um      non_rigid_D 20.50 -> 3.855 um
    6442   rigid_D  6.89 -> 3.331 um      non_rigid_D 16.33 -> 3.414 um

Note the second column. The non-rigid stage used to make VALIS's own metric 1.3-2.4x **worse**
than the rigid fit it started from; it no longer does.

**2. The micro-registration pass was inert.** `register_micro` was called without
`reference_img_f`, so it forwarded `None` to `serial_non_rigid.register_images`, where
`warp_tools.get_ref_img_idx` returns index 0 for two images (warp_tools.py:57-59) — the **Xenium**.
That made the Xenium the fixed image, so it was handed a field of zeros
(serial_non_rigid.py:778-781) while the field computed for the IF was skipped by the update loop,
which skips the *Valis object's* own reference — the IF (registration.py:4941-4943). The two
references disagreed, the computed field was discarded, and the only surviving effect was that the
**old** field got resampled: measured, `non_rigid_D` went 20.497 → 21.331 on 6436, i.e. slightly
worse. Passing the argument makes the pass do what its name says. It would not have fixed the
residual on its own in any case: every path in `register_micro` ends at
`bk_dxdy = current + new` and none of them writes `M`, so it cannot remove a translation.

**3. Neither pyramid is `micron = pixel * um_per_px`, and they are wrong differently.** Measured on
donor 6436, three 384-px tiles, level 0 against level 2:

| | box mean | best decimation |
|---|---:|---:|
| Xenium `morphology_focus` | 0.9766 | **0.9936** at phase (0,0) |
| IF qptiff (DAPI page) | **0.99992** | 0.9891 |

So the Xenium pyramid *keeps* level-0 pixel `f·k` and the IF pyramid *averages* `f·k … f·k+f−1`.
Their pixel centres sit at `(f·k + 0.5)·base` and `(f·k + f/2)·base`, i.e. `origin_um` of 0.10625
and 0.5. Ignoring both left a fixed **+0.39 µm per axis** in the warped result. Separately,
`shapes[0][0] / shapes[lvl][0]` is not the downsample factor: on 6436 it is 68119/17029 = 4.000176
where the factor is exactly 4, a scale error worth 0.61 µm at the far corner. Both are now derived
from the integer factor with an assertion, and both are recorded in `export.json`.

## Every candidate, scored the same way

All of them start from the same seed and are scored by the same per-cell distance, and the losers
are reported rather than dropped (`refine_candidates.png`):

| candidate | package | 6436 | 6442 |
|---|---|---:|---:|
| `valis_affine` — the seed | valis-wsi 1.2.0 | 1.190 | 1.830 |
| `valis_field` — affine + its 250 µm binned residual | valis-wsi | 1.486 | 1.889 |
| `nuclei_translation` — seed + one constant | scipy | 0.801 | 1.809 |
| `nuclei_icp_similarity` | `register.icp` | 0.791 | 1.805 |
| **`nuclei_icp_affine`** | `register.icp` | **0.546** | **0.630** |
| `image_translation` — whole-section NCC | numpy FFT | 0.801 | 1.809 |
| `phase_translation` — phase correlation, upsample 100 | scikit-image | 2.259 | 1.891 |
| `nuclei_icp_affine+local` — plus a 250 µm field | scipy | 0.522 | 0.600 |

Four things this table says that a single number would not:

* **The affine matters and is not shear-hunting.** It cuts the residual 1.4x on 6436 and 2.9x on
  6442 against the similarity fit, and it does so by *reducing* anisotropy — 643 → 280 ppm on
  6436 and 1338 → 419 ppm on 6442. It is removing a stretch VALIS had, not adding one.
* **The local field does not earn its place.** It buys 0.02-0.03 µm, inside the declared 0.05 µm
  tolerance, so `choose()` keeps the simpler transform. More freedom cannot score worse on the
  data it was fitted to, which is exactly why the tolerance is declared rather than "lowest wins".
* **`valis_field` is now worse than `valis_affine` on both donors.** A 250 µm median-binned
  residual field with a third of its bins filled from a neighbour was compensating for a bad
  affine; against a good one it only adds noise.
* **Three mechanisms were asked for the same translation.** Nearest-neighbour nuclei and
  whole-section NCC agree to **0.046 µm (6436)** and **0.092 µm (6442)** — two different
  mechanisms, one answer. Phase correlation is the outlier, 2.15 and 0.53 µm from the other two,
  and it is reported rather than hidden: normalising the cross-power spectrum weights every
  spatial frequency equally, and a sparse point raster against a continuous DAPI image is exactly
  the case where that hurts.

SimpleITK 2.5.6 (in `~/.venvs/valis`) and itk-elastix 0.19.2 (in `~/miniconda3/envs/image2image`)
are installed and were **not** run. They would refine the same seed by the same intensity mechanism
as the two image candidates above, which the point-set fit beats by 1.5x (6436) and 2.9x (6442)
on the same score.
That is a judgement, and it is recorded here so it can be reversed.

## Is it a registration, or a fit?

Mutual nearest-neighbour matching can be self-consistently **wrong**: if the seed were out by more
than half the spacing between nuclei, every point would pair with a stranger and the fit would
converge confidently onto that, with nothing in the in-sample score showing it. So the section is
split in two and each fit is scored on the half it never saw:

| donor | fitted on | in-sample | **held out** |
|---|---|---:|---:|
| 6436 | left half (348,232 nuclei) | 0.542 | **0.614** |
| 6436 | right half | 0.545 | **0.585** |
| 6442 | left half (479,459 nuclei) | 0.657 | **0.605** |
| 6442 | right half | 0.603 | **0.671** |

Held-out is within 0.07 µm of in-sample, and on 6442 one fold's held-out score is *better* than its
in-sample one. The two halves' transforms disagree by 0.184 µm (6436) and 0.240 µm (6442) of
translation and by 19 and 2 ppm of scale. The transform transfers, so it is a property of the
sections rather than of the points it was fitted to.

## The gate, and what it is now worth

`run_align.py` warps Xenium nucleus **polygons** into the IF frame, rasterises them at 0.5 µm/px and
finds, per 500 µm tissue window, the shift that best lines that mask up with the IF DAPI by
normalised cross-correlation. Its thresholds are geometry, not preference: PASS at the displacement
where a disc of the donor's own median cell still covers 75 % of its own footprint.

| donor | windows | median | p90 | precision | DAPI contrast inside vs gaps | verdict |
|---|---|---:|---:|---:|---:|---|
| 6436 | 243/243 | **0.109 µm** | 0.233 | 0.015 | 2.04 → **5.14×** | PASS (bar 1.795) |
| 6442 | 365/365 | **0.134 µm** | 0.260 | 0.014 | 3.51 → **4.98×** | PASS (bar 1.880) |

**Read this as consistency, not independence.** For the tissue-outline and VALIS transforms the
gate was genuinely independent — neither had ever seen a nucleus. The nucleus-matched fit and the
gate now optimise closely related things, so the gate confirming at a tenth of a micron is expected.
The checks that are *not* circular are the hold-out above and the placement null.

**Which diameter anchors the bar is a decision, and all three are reported.** The gate takes its
diameter from `cell_area`, but it correlates *nucleus* polygons, and the hormone lane samples an
8 µm disc:

| | 6436 | 6442 |
|---|---:|---:|
| cell (9.08 / 9.51 µm) — the declared bar | 1.795 | 1.880 |
| nucleus (6.19 / 6.13 µm) — what the gate actually matches | 1.223 | 1.211 |
| the hormone sampler's 8 µm aperture | 1.581 | 1.581 |

The per-cell median (0.546 / 0.630 µm) is under all three, and so is the per-cell p90 for the cell
and aperture bars. The declared cell bar is what `run_align.py` reports, because it was declared
before the run and is not renegotiated in either direction.

**Tissue Dice went slightly down**: 0.8903 → 0.8813 (6436) and 0.9509 → 0.9497 (6442). That is not a
regression, it is the point. The ECC transform was *fitted to maximise* tissue-outline overlap, so
it wins on Dice by construction while sitting 42 µm out at cell scale. Dice is not the criterion.

## Where the remaining 0.55 µm is, and where it is not

* **Not a translation.** After the affine the median residual vector is (−0.006, +0.000) µm.
* **Not spatial structure.** A 250 µm field with full freedom to bend changes the median by 0.02 µm,
  and the residual quiver has no direction anywhere on the section (`refine_residual.png`).
* **Not peak localisation.** Detecting the IF nuclei at 0.5 µm/px instead of 1.0 µm/px moves the
  residual to 0.563 (6436) and 0.650 µm (6442) — very slightly *worse*, not better. Doubling the
  resolution buys nothing, so the limit is not how precisely a peak can be located.

What is left is the difference between where Xenium's nucleus-polygon centroid sits and where the
IF DAPI blob centre sits: two definitions of one nucleus, plus whatever genuinely differs between
two imagings of the same section. No transform of any class removes that.

**0.5 % of Xenium nuclei have no IF nucleus anywhere near them.** They are outside the drawn IF
`Tissue` polygon — on 6436 the drawn tissue stops 1.3 mm short of where the Xenium cells reach —
so there is nothing there to register to. Those are recorded as far distances rather than dropped,
and they are what the tail of `refine_pairs.png` is.

## The `Islet` class needs a size filter, and it is not optional

The classifier emits thousands of sub-cell fragments beside the real islets:

| donor | objects | median area | ≥50 µm across | share of islet area |
|---|---:|---:|---:|---:|
| 6436 | 4,671 | 27 µm² (5.9 µm) | **125** | 84.6 % |
| 6442 | 3,366 | 23 µm² (5.4 µm) | **95** | 81.1 % |

The median object is one cell wide. Filtered, the counts land in the range the PhenoCycler cohort
shows (67–1,304 islets per donor).

## The hormone check, re-run on this transform — and it still does not reach the bar

`run_hormone.py --transform refined` → `hormone_panc_refined.json`. Same AUROC, same null, same
pre-registered bar of 0.75; only the registration changed. (The `ecc` numbers on disk carried two
defects at once — 31–42 µm of misregistration *and* `um_per_px = 1.0018796`, the qptiff resolution
tag rather than the recorded 0.25 — so they were never a usable "before". The `base_um_per_px`
guard now in `read_qptiff_channel` closes the second one.)

| | INS ranks beta-like | | GCG ranks alpha-like | |
|---|---:|---:|---:|---:|
| | ecc (31–42 µm) | **refined (0.55 µm)** | ecc | **refined** |
| 6436 | 0.533 | **0.514** | 0.517 | **0.522** |
| 6442 | 0.574 | **0.657** | 0.562 | **0.728** |

**The bar is not met on either donor and it is not moved.** What the two rows say is not one
result but two, and the difference is the point:

* **6442 moved a long way and in the right direction** — GCG 0.562 → 0.728, and
  `agreement_by_scale` at zero radius went **+0.074 → +0.486**. A registration fix producing that
  is a registration that was previously costing real signal.
* **6436 did not move at all** (0.533 → 0.514), and its zero-radius agreement is +0.055. Its
  registration is not the reason: `bundle_qc.png` shows its alignment peaking at zero offset
  exactly as 6442's does.

**The first reading of that — "6436's staining does not distinguish islets" — was wrong, and
`channel_check.py` is what corrects it.** 6436 has an excellent islet channel. It is **Opal 520**,
the one recorded GCG, at **66× inside islets against elsewhere**; and it is **Opal 690**, the one
recorded INS, that carries no islet signal (2.2×, p99.9 of 15 counts against 98 for Opal 520).
Measured in the qptiff's own frame against the maintainer's own polygons, with no registration in
the path at all:

| donor | DAPI | Opal 570 (FAILED) | Opal 690 (INS) | Opal 520 (GCG) |
|---|---:|---:|---:|---:|
| 6436 | 1.4× | 1.9× | **2.2×** | **66.0×** |
| 6442 | 1.1× | 117× | **17.1×** | **401×** |

On 6442 both hormone channels are islet-restricted and visibly complementary *within* one islet,
which is what INS and GCG should look like. On 6436 exactly one of the two worked, and it is not
the one the AUROC was reading as INS. Whether that means the recorded mapping is wrong for this
section or that its INS stain failed is **the maintainer's recorded fact to resolve, not a
measurement** — `restain_sections.csv` is hand-recorded precisely because it is in no file.

Note also, from the same sweep across all five recorded sections: `Opal 570` is recorded FAILED on
every row, but on donors 115 (0041117) and 6591 (0041116) — the other acquisition batch — it
carries a p99.9 of 253 and 156 against 0–1 on the three nPOD-batch sections. The blanket FAILED
recording is a nPOD-batch fact that may not extend to those two.

`agreement_by_scale` still climbs monotonically to 120 µm on both donors, which remains the
reason to read this lane as supporting islet-level cross-modality claims and not per-cell ones.

## `validate.py` has NOT been re-run

*IF islets against Xenium endocrine RNA* (`validate.py`, `registration_panc.json`): 6442 17.30×
inside vs elsewhere, 6436 2.83×, against a placement null of ~0.117. 6436's alpha and beta moved
together (2.83× / 2.89×) where 6442's did not (6.33× / 25.38×), and `validate.py`'s own docstring
calls that the signature of a mask landing partly off-target. Those numbers were measured on the
**tissue-outline** transform and **must not be read as current**.

## Historical biological Explorer bundle (6436 and 6442)

This section describes the older, antibody-role-dependent `run_explorer_export.py` delivery for
the two reviewed donors. It is retained as provenance for those protected bundles; it is **not**
the cohort registration-only export described above. `run_explorer_export.py` writes, per donor, a
directory that needs neither this repository nor this environment:
`data/development/restain/explorer_bundle/donor_<id>_<serial>/`.

| file | what it is |
|---|---|
| `images/<serial>_IF.ome.tif` | DAPI + INS + GCG, pyramidal, tiled 1024, 0.25 µm/px, lossless — 10x's stated requirement for an importable post-Xenium image. 2.1–2.5 GB. Opal 570 (FAILED) and Sample AF are **not in the file**; `manifest.json` records that they exist and why they were dropped |
| `images/<serial>_IF_imagealignment.csv` | the 3x3 row-major matrix, no header, last row `0,0,1`, mapping **IF image pixels → Xenium morphology pixels**, in the byte layout of Explorer's own export (below). Also valid input to `xeniumranger import-segmentation` and `spatialdata_io.xenium_aligned_image` |
| `images/<serial>_IF_similarity_imagealignment.csv` | a **diagnostic**: the nearest similarity to the affine, the only shape Explorer's own aligner ever exports; its cost against the affine is in `manifest.json → similarity_diagnostic` |
| `images/<serial>_IF_keypoints.csv` | `fixedX,fixedY,alignmentX,alignmentY`: 20 Xenium nuclei and their exact pre-images under the matrix — Explorer's other documented input, which it fits its own similarity through |
| `images/<serial>_IF_on_xenium_grid.ome.tif` | the same IF resampled onto the Xenium morphology grid (same width × height, 0.2125 µm/px, same 9 levels) through the same transform, for **"Yes, skip alignment"** — no CSV involved |
| `cell_groups/hormone_positivity.csv` | `cell_id,group,color` for all 695,113 (6436) and 958,744 (6442) cells: INS/GCG positivity at the recorded threshold |
| `README.md`, `manifest.json` | the instructions (which file goes with which Explorer option, in which order, and what to report) and the provenance; `bundle_qc.png` is regenerated beside them by `explorer_figures.py` |

**The direction of the matrix is the one thing here that cannot be recovered by a careful
reader**, so it is derived rather than assumed. This lane's affine maps Xenium µm → IF µm;
Explorer wants IF px → Xenium px, so the matrix is `S_xen⁻¹ · A⁻¹ · S_if`. Three sources agree on
the target: 10x's file format doc (0.2125 µm/px, origin at the image's upper left, "divide by the
pixel size"), 10x's alignment doc ("row-major", "units are in pixels"), and `spatialdata_io`,
which reads the CSV `header=None` and applies it as an `Affine` into the frame it puts Xenium
geometry in by scaling microns by `1/pixel_size`. `S_if` carries the **measured** box-mean pixel
origin; dropping it is a flat 0.125 µm per axis, and a test plants exactly that.

One convention is genuinely unresolved and is recorded rather than assumed away: no 10x document
says whether Explorer indexes an imported image by pixel corner or centre. That is ±0.125 µm.

**Why there are four ways in, and what each one is for.** Explorer 4.x refused the first
version of the matrix file with "Failed to load file." and showed the native image unaligned
under "Yes, skip alignment". The second is not a defect: 10x documents that option only for an
image already resampled into the Xenium coordinate system (its reference file for it has
exactly the morphology image's dimensions at 0.2125 µm/px), and Explorer places such an image
with the identity in pixel space, ignoring OME PhysicalSize. The first has no documented cause —
10x publishes no validation rule for the import — so everything in which the file differed from
Explorer's own export was removed, and the other documented inputs were added:

* **Bytes.** Three 10x-published `_imagealignment.csv` files were hexdumped. Explorer writes
  CRLF between rows, **no newline after the last row**, the last row literally `0,0,1`, and at
  most nine significant digits. The old writer (`csv.writer` + `repr`) added a trailing CRLF
  and wrote `0.0,0.0,1.0` with 17-digit floats. The new writer reproduces the 10x exemplar
  byte for byte (a test feeds it the exemplar's numbers and demands the same 73 bytes); nine
  digits move the far corner of the field by under 3e-4 Xenium px.
* **Structure.** Every Explorer-exported matrix is an *exact* similarity (`a == d`, `b == -c`);
  Explorer's aligner is documented as similarity-only. Ours is a general affine (anisotropy
  2.8e-4 / 4.2e-4), and whether the importer minds is unknown. So the bundle also carries the
  nearest similarity (through the keypoints, canonicalised so the equalities hold bit-exactly)
  and the keypoints themselves, which let Explorer fit what it accepts and keep the result
  adjustable. Neither is a re-fit: both are the shipped matrix re-expressed at chosen points,
  and the manifest records what they cost (`similarity_diagnostic.deviation_from_affine_um`).
* **The warped image.** `<serial>_IF_on_xenium_grid.ome.tif` is the IF resampled onto the
  morphology grid: output pixel `(k, l)` samples the native image at `M⁻¹(k + ½, l + ½)`
  (Xenium pixels are corner-origin; the matrix's IF side is centre-indexed; cv2 treats integer
  coordinates as centres — all three pinned by tests), bilinear on the stored 8-bit values,
  zeros outside the IF, and a floor-halved box-mean pyramid with exactly the morphology's level
  shapes. It is the one delivery whose placement depends on nothing but documented Explorer
  behaviour. Positivity is still measured on the native raster, never on this one.

Ruled out along the way: CRLF itself and float precision as *definite* causes (Explorer's own
files use CRLF, and two other 10x files carry 17-digit floats — whose producer is unknown, so
precision stays matched rather than argued), image size (10x ships Explorer-aligned images
larger than ours on every metric), missing TIFF resolution tags and 1024-px tiles (10x's own
skip-alignment file has neither), matrix direction, and pyramid level.

**`bundle_qc.png` checks the bundle the way the colleague will read it** — from the shipped
files only (native OME-TIFF, matrix CSV, similarity CSV, keypoints CSV, warped OME-TIFF), never
from the transform that made them, so it is able to disagree with them. Three routes into the
same overlay, each read as Explorer would read it: the native image through the shipped matrix,
through the shipped similarity, and the warped image with the identity (`px = µm / 0.2125`,
floor). Mean DAPI at Xenium nucleus centroids against a planted offset:

| donor | route | at 0 µm | between nuclei | ratio | argmax |
|---|---|---:|---:|---:|---:|
| 6436 | native + matrix | 49.52 | 22.38 | 2.21× | **0 µm** |
| 6436 | native + similarity | 49.45 | 22.37 | 2.21× | **0 µm** |
| 6436 | warped + identity | 49.49 | 22.38 | 2.21× | **0 µm** |
| 6442 | native + matrix | 39.76 | 19.08 | 2.08× | **0 µm** |
| 6442 | native + similarity | 39.71 | 19.09 | 2.08× | **0 µm** |
| 6442 | warped + identity | 39.72 | 19.08 | 2.08× | **0 µm** |

Both curves fall to the floor by 4 µm and dip *below* it at 6 µm — one nucleus radius out lands in
the gap between nuclei. A merely-plausible alignment does not produce that shape, and the warped
image read with the identity reproduces the native curve to the second decimal (peak ratio
0.999, 0.999).
Warp cross-check on the shipped files: 6436: p99 |diff| 0.00, max 1 count, 100.0 % within 1 count; 14/14 checks, PASS; 6442: p99 |diff| 0.00, max 1 count, 100.0 % within 1 count; 14/14 checks, PASS.

The warp is pinned a second way: a 512² Xenium-px window of the shipped warped DAPI is
differenced against the native DAPI resampled with `scipy.ndimage.map_coordinates` at
`M⁻¹(k + ½, l + ½)` — an implementation sharing no code with the writer. A half-pixel error
would show as tens of counts at nucleus edges; the writer's 1/32-px fixed point leaves about
one. `bundle_qc.json` carries every number and a PASS/FAIL per donor.

## What the bundle refuses to claim

*Stale as of 2026-09-02: this section describes the deciles-only bundle of `9e88272`. The current
driver ships one absolute positivity cut (20 of 255, mean inside the cell's own polygon), recorded
in the manifest and changeable with `--threshold`; it is left here unrewritten because the choice
between the two is the maintainer's.*

The maintainer's rule for this export is that a per-cell hormone-positive label ships only if the
check it rests on clears the bar that check was pre-registered at. Neither marker does, on either
donor, so **no cell in either bundle is called INS+ or GCG+**. Both channels ship as deciles of
the measured intensity, per donor, on that donor's own distribution, and the bundle README says in
its own words that a `d10` cell is a cell where a channel read high and nothing more.

The other branch is implemented and tested, not hypothetical: `_hormone_groups(..., cleared=True)`
fits a two-component mixture to the donor's own log intensity and emits `INS_positive` /
`INS_negative`. Which branch runs is read out of `hormone_panc_refined.json`. A test plants the
same bimodal data at mixing weights of 5 %, 25 % and 50 % and requires the threshold to move by
less than 0.5 in log space — the fit is not steered by how many cells come out positive, which is
what the repository's no-frequency rule actually requires.

## Files

| file | role |
|---|---|
| `registration_inputs.py` | strict registration inventory loader, qptiff facts, Xenium cross-checks, and central donor eligibility |
| `export_dapi.py` | exact IF/Xenium DAPI planes, pyramid factors, pixel origins, and input provenance |
| `nuclei.py` | annotation-optional IF nucleus detection and Xenium nucleus centroids from the bundle |
| `refine.py` | the candidates, the per-cell score and its placement null |
| `run_refine.py` | manifest-backed driver, parsimony rule, spatial hold-out, donor-local `refinement.json`, and `refined` transform |
| `refine_figures.py` | candidate comparison, per-cell distribution, residual quiver |
| `valis_register.py` | the seed (VALIS venv) → `warped_if_um.npy` |
| `valis_transform.py` | VALIS's point map as a `Transform`, for the comparison candidate |
| `align.py` / `run_align.py` | all-window or held-out-fold DAPI gate → donor-local `alignment_<stem>.json` |
| `run_intensity_refine.py` | optional one-fold robust DAPI correction → new exact transform stem and fit report |
| `align_figures.py` | the gate's own QC sheets |
| `run_registered_export.py` | PASS-only, exact-stem, transactional export of every raw channel onto the Xenium grid; no biological calls |
| `islet_annotations.py` | strict QuPath-pixel → IF-micron → full inverse transform → Xenium-grid-pixel annotation writer and Explorer validation |
| `run_islet_annotation_export.py` | annotation-only transactional publisher for an already-rendered Explorer bundle, including protected reviewed images |
| `pheno_annotations.py` | serial-section PhenoCycler-pixel → micron → S4 inverse → Xenium-grid-pixel multi-class writer, field clipping, and over-cap splitting |
| `run_pheno_annotation_export.py` | annotation-only additive publisher for the PhenoCycler classes, transactional against an existing bundle |
| `pheno_annotation_figures.py` | overlay QC sheet, drawn only from the shipped sidecars and the Xenium cell table |
| `pheno_rotation_figures.py` | rotation-response diagnostic: out-of-sample islet Dice/RMSE beside footprint coverage |
| `geo.py` | the exports in microns, the islet size filter, the scale guard |
| `tiles.py` / `rasterize_labels.py` / `display.py` | windowed reads, rasterisation, honest ranges |
| `validate.py` / `hormone.py` / `run_hormone.py` | historical biological downstream checks (not part of registration-only cohort export) |
| `run_register.py` | historical annotation-derived tissue-outline registration; also scores a supplied transform |
| `explorer.py` | shared morphology-grid, pyramid, OME-TIFF, warp, matrix, keypoint, and legacy cell-group primitives |
| `run_explorer_export.py` | historical antibody-role-dependent bundle driver for the two protected reviewed donors |
| `explorer_figures.py` | historical biological bundle QC, drawn from shipped files |
| `tests/` | synthetic geometry, inventory, fold-isolation, exact-transform, and registered-OME contract tests |

`data/development/restain/valis_baseline/` holds the 2026-08-31 arm this replaced.
