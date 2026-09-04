# DAPI / Phenocycler Follicle U-Net

> **Input is the DAPI channel alone.** Follicle detection is a nuclear-density
> problem, so DAPI carries most of the signal that the deconvolved hematoxylin
> channel carries in the H&E pipeline — which is exactly why the two cohorts can
> be pooled. What DAPI does *not* carry is the marker information the training
> labels were built from; see [The DAPI-only ceiling](#the-dapi-only-ceiling).

Companion to `HE_unet_README.md`. Same task, same 2 µm/px grid, same donor-level
CV, same success bar. Read that one first if you have not.

---

## Why this exists

The H&E cohort's follicles were hand-drawn; the U-Net replaces that manual step.
The fluorescent cohort is the mirror image — its follicles were already
segmented, by the `SpleenRegions*` QuPath pixel classifier in
`scripts/Workflow.groovy`, using 25–33 marker channels. That gives us **14 more
slides of free supervision** for the same task.

So this pipeline does two things:

1. **Distils** the multi-channel pixel classifier into a DAPI-only model. The
   student sees one channel; the teacher saw thirty. Whatever the student
   recovers is a statement about how much of the region structure is visible in
   nuclear density alone.
2. **Feeds the H&E model.** Both pipelines write single-channel images where
   *high nuclear density is dark*, in the same directory layout, so
   `train_HE_unet.py --extra-data-dir data_dapi` pools them into one training set
   without any conversion step.

---

## Files

| File | Role |
|---|---|
| `export_DAPI_training_pairs.groovy` | QuPath → `data_dapi/{images,tissue,masks,regions,meta}` at 2 µm/px |
| `DAPI_unet_utils.py` | Normalisation parity fns, multiplex reader, multitask dataset/loss/metrics |
| `train_DAPI_unet.py` | Donor-level k-fold training → `analysis/DAPI_unet_models/` |
| `predict_DAPI_unet.py` | Whole-slide inference → GeoJSON + `Measurements/DAPI_follicles_unet.csv` |
| `check_dapi_parity.py` | Verifies the export and inference transforms agree **to the byte** |
| `report_DAPI_concordance.py` | Per-donor concordance vs the pixel classifier, incl. Follicle↔PALS confusion |
| `make_synthetic_dapi_data.py` | Synthetic donors + synthetic multiplex OME-TIFFs for tests/CI |
| `../tests/test_DAPI_unet_utils.py` | Unit tests, including Groovy-transcription parity checks |

---

## Quickstart

```bash
# 1. Deps (same file as the H&E pipeline — no new dependencies)
pip install -r scripts/requirements_HE_unet.txt

# 2. Export training data from the QuPath project (read-only; no --save needed)
QUPATH="$HOME/QuPath-v0.7.0-Linux/QuPath/bin/QuPath"
"$QUPATH" script --project project.qpproj scripts/export_DAPI_training_pairs.groovy

# 3. Verify the export and the inference path agree, BEFORE training on it
python scripts/check_dapi_parity.py --data-dir data_dapi --images-dir Images

# 4a. Find the plateau first (single fold, generous patience)
python scripts/train_DAPI_unet.py --folds 1 --epochs 60 --patience 12

# 4b. The reportable run
python scripts/train_DAPI_unet.py --folds 4 --epochs 50

# 5. Honest reporting
python scripts/report_DAPI_concordance.py

# 6. Whole-slide inference
python scripts/predict_DAPI_unet.py --wsi-dir Images \
    --models 'analysis/DAPI_unet_models/fold_*.pt'
```

Optional exporter arguments:

```bash
"$QUPATH" script --project project.qpproj \
    -a TARGET_UM_PER_PX=2.0 -a EXCLUDE=HDL018,HDL021 -a OUT_DIR=data_dapi \
    scripts/export_DAPI_training_pairs.groovy
```

---

## Design choices

### 2 µm/px — set by the labels, not by taste

`SpleenRegions2.json` declares `inputResolution = 4.062131048194629 µm`, an exact
**8× downsample** of the 0.5078 µm/px Phenocycler grid (`SpleenRegionsCODEX` is
2.992 µm on the 0.374 µm grid — also 8×). The region boundaries therefore carry
**no detail finer than ~8 µm** and were upsampled to native coordinates.

Training at native 0.5 µm/px would mean learning 8×-blocky targets at 16× the
cost. 2 µm/px is still twice as fine as the labels, and it is the grid the H&E
pipeline already uses — which is what makes pooling possible.

`SmallVessel` is the opposite case: labelled at native 0.508 µm/px with a 10 µm²
floor (~3.5 µm across), i.e. sub-pixel at 2 µm/px. It is excluded from the region
classes. A vessel model needs its own resolution and its own script.

### Single-channel, inverted, per-slide normalised

Nothing upstream normalises these images — the pixel classifiers split on raw
uint16 (thresholds like `le: 21871.7`), so their labels carry per-slide staining
bias that the input should not inherit. The export therefore normalises per
slide, and **inverts**:

```
bin      = value >>> 8                                   # 256-bin histogram
t1, t2   = 3-class Otsu                                  # t1 = glass/tissue
tissue   = raw >= ((t1 + 1) << 8)
lo, hi   = 1.0 / 99.5 percentile of raw WITHIN tissue
intensity= floor(255 * (1 - clip((v - lo)/(hi - lo), 0, 1)))
```

Two details worth their own line:

* **Why three-class Otsu.** A spleen DAPI histogram has *three* modes — glass,
  red pulp, white pulp. Two-class Otsu puts its threshold between red pulp and
  white pulp, producing a "tissue" mask containing only white pulp, which would
  then supply every negative training patch and every percentile. The 3-class
  lower threshold is the glass/tissue split. The sidecar records `tissue_frac`
  and the export warns below 5% so a bad split is visible immediately.
* **Why inverted.** `HE_unet_utils.rgb_to_h` writes high hematoxylin as *dark*.
  Raw DAPI is high nuclei = *bright*. Pooling two cohorts that disagree about
  which way "more nuclei" points would be worse than not pooling at all.

### The parity contract

The transform above is implemented **twice** — in Groovy for the export, in
Python for inference. If they drift, the model trains on one intensity scale and
predicts on another, and nothing else in the pipeline notices.

Every step is integer or `double` arithmetic over a 256-bin histogram
specifically so the two languages agree to the byte. That is defended at two
levels:

* `tests/test_DAPI_unet_utils.py` contains **imperative transcriptions** of the
  Groovy loops (`_reference_otsu`, `_reference_otsu3`, `_reference_percentile`,
  `_reference_normalise`) and asserts the shipped vectorised code matches them.
  Written as loops on purpose — agreement is evidence, not a tautology.
* `check_dapi_parity.py` re-derives every exported PNG from its OME-TIFF through
  the Python path and compares. It needs no QuPath, so it runs anywhere the
  images live, and CI runs it against synthetic OME-TIFFs.

The H&E pipeline documents its equivalent check as a manual README snippet; this
one is automated. **Any change to the constants or the arithmetic must be made in
both files in the same commit.**

**Why the pixel comparison has a tolerance while the parameters do not.** The
two sides reach the 2 µm/px grid by different routes: QuPath resamples with its
own filter at a possibly non-integer downsample (3.94 for the 0.5078 µm/px
slides), while Python reads the nearest suitable pyramid level and resizes with
PIL. Both are correct; they are not bit-identical, so raw counts can differ
slightly. After normalisation that is invisible — a ~50 000-count window maps to
256 levels, so ~200 raw counts per output level — but not guaranteed to be zero.

So `otsu_bin`/`lo`/`hi` are compared **exactly** (integer bin indices and bin low
edges, unmoved by resampling noise — a difference there is genuine drift), and
the images are compared with `--max-abs-diff` / `--max-mismatch-frac`. On
synthetic data both paths read the same raster with no resampling between them,
so CI runs `--strict` and demands byte-exactness. The transform arithmetic
itself is pinned exactly, with no resampling in the way, by the unit tests.

### Multitask head

One decoder, `classes = 1 + K`:

* **channel 0** — binary Follicle logit, `0.5·BCE + 0.5·Dice`. This is the
  deliverable and the *only* thing model selection is scored on, so a checkpoint
  can never win on region accuracy while the follicle segmentation degrades.
* **channels 1..K** — region softmax over
  `Background / RedPulp / Trabeculae / LargeVessel / PALS / Follicle`,
  weighted CE + macro-dice, at `λ = 0.3`.

The auxiliary loss is **per-sample masked**: donors with no `regions/` map
contribute zero auxiliary gradient. That is the entire pooling mechanism — an
H&E donor trains the Follicle head and is silently invisible to the region head.

Two details that matter for the numbers:

* Class weights are `1/sqrt(frequency)` over the **training** donors, computed
  once and stored in the checkpoint. Unweighted, the region head would score well
  by predicting RedPulp nearly everywhere. Classes absent from the training set
  get weight **0**, not a huge one.
* Macro-dice averages over classes **present in that patch's target**. A 512 px
  patch rarely contains all six regions, and charging a full 1.0 dice penalty for
  the absent ones makes the loss depend mostly on what happened to be in frame.

`--no-aux` gives `K = 0`, a model architecturally identical to the H&E one (there
is a test asserting the state dicts have identical shapes), so checkpoints
transfer in either direction.

### Region label space and paint order

```
0 Background/Unlabeled   1 RedPulp   2 Trabeculae   3 LargeVessel   4 PALS   5 Follicle
```

The exporter paints in that order, later winning, so **Follicle is painted last**
and the binary `masks/` PNG is *derived from* the label map — the two cannot
disagree. `check_dapi_parity.py` re-asserts it, and so does the test suite.

`Ignore*` produces no annotation (`Workflow.groovy` does not pass
`INCLUDE_IGNORED`), so index 0 means "the classifier declined to label this", not
"confidently background". Worth remembering when reading the aux head's numbers.

One Java2D trap, since it is easy to reintroduce: the label map is painted into
`TYPE_INT_RGB`, **not** `TYPE_BYTE_GRAY`. Filling a grey image with
`new Color(3,3,3)` runs an sRGB→linear-grey conversion and does *not* land on
sample value 3. The H&E exporter only uses black and white, so it never hit this.

### Donor-level splits, and a leak that was already there

`donor_kfold` groups **files by donor before dealing**, using the explicit
`sample_id` loaded from each metadata JSON or resolved by exact lookup in
`config/sample_manifest.csv`. This fixes two historical failure modes:

* Multiple exports declared as `HDL073` remain in one fold rather than being
  dealt round-robin as separate files.
* `1901HBMP004` is genotype-stratified by its declared identity, without
  interpreting any filename prefix.

14 entries → 14 donors. `--folds 4` leaves ≥3 validation donors per fold, which
is what the pooled Spearman needs; `--folds 5` leaves 2 in some folds and the
trainer prints a note when that happens.

---

## The DAPI-only ceiling

**Follicle vs PALS is not decidable from DAPI alone.** The teacher separated them
with CD20 vs CD3e. Both are dense-nuclei white pulp. A DAPI-only model can post a
respectable follicle *count* while systematically swapping the two zones, and a
report led by mean IoU would hide exactly that.

So the reporting is deliberately ordered:

1. per-class region IoU
2. the **Follicle→PALS / PALS→Follicle** confusion rates
3. only then the headline count/area Spearman

`--merge-white-pulp` collapses the two into one `WhitePulp` class. If scores jump,
the residual error really was that boundary rather than segmentation quality —
run it as an ablation, not as the default.

A class absent from both prediction and truth reports `None`, **not** 1.0. Several
donors have no Trabeculae or LargeVessel at all, and scoring those as perfect
would inflate the mean with classes the model was never asked to find. (The
binary `iou_score` keeps its both-empty == 1.0 convention; that is a different
question.)

### Success bar

Same as H&E, and it is **not IoU**: per-donor follicle **count + mean-area
Spearman ρ > 0.9** against the teacher, from the `--folds 4` run
(`pooled_metrics.json`). Model selection uses semantic IoU, so a high IoU can
still mis-count after watershed — which is why instance metrics are reported
separately.

### Train to the plateau

The H&E model's previously reported IoU ≈ 0.715 was an artifact of `--epochs 5`
making the best epoch coincide with the cosine LR annealed to ~1e-6 (recorded in
`CLAUDE.md`). The same trap applies here. Run `--folds 1 --epochs 60
--patience 12` first and look at the per-epoch history that `cv_summary.json`
embeds before committing to a reportable run.

---

## Pooling the two cohorts

Because both pipelines write the same contract, pooling is a flag:

```bash
# DAPI model, pooled with H&E donors (aux head sees DAPI donors only)
python scripts/train_DAPI_unet.py --data-dir data_dapi --extra-data-dir data

# H&E model, pooled with DAPI donors
python scripts/train_HE_unet.py --data-dir data --extra-data-dir data_dapi

# H&E model warm-started from a DAPI checkpoint (head shapes differ; only
# matching tensors load, and the count of skipped ones is printed)
python scripts/train_HE_unet.py --init-from analysis/DAPI_unet_models/fold_0.pt
```

Donor names must be unique across pooled dirs — a collision aborts rather than
being renamed, because it would corrupt the dataset's name-keyed cache *and* hide
a train/val leak at the same time.

**The ablation to actually run** is H&E-only vs pooled vs DAPI-pretrained, scored
on the same H&E folds. Pooling is a hypothesis here, not a result: it assumes
inverted-normalised DAPI and deconvolved hematoxylin are close enough in
appearance for shared weights to help. Measure it before believing it.

---

## Whole-slide inference

`read_dapi_channel` goes through **tifffile and the pyramid levels**, not
OpenSlide — OpenSlide cannot select one channel of a 33-channel uint16 OME-TIFF.
It resolves the DAPI channel by OME-XML name and falls back to index 0 with a
warning (channel names have been lost on some of these images before, which is
why `RestoreChannelNames.groovy` exists), then picks the finest level that is
still at least as coarse as the target and resamples to exactly 2 µm/px.

Only the **Follicle channel** is stitched across a slide. The model still
computes every channel; `follicle_activation` just stops `predict_wsi` from
accumulating the rest. A slide at 2 µm/px is ~10⁸ px, so stitching 1+6 float32
channels costs ~3 GB against ~430 MB for the one channel anything downstream
reads — and the same saving applies to the per-epoch validation, which runs on
full-slide rasters every epoch.

Region metrics do need all K channels, so the trainer computes them
**band-wise** (`region_confusion_banded`): each band owns a disjoint range of
rows, is predicted with a one-patch margin for context, and contributes to a
running confusion matrix. Confusion matrices add, so the result matches a single
pass while peak memory stays around 1 GB.

**GeoJSON units.** `predict_DAPI_unet.py` defaults to `--geojson-units px` so the
output re-imports into QuPath via *Objects → Import objects*. The existing
`analysis/HE_unet_predictions/*.geojson` are in µm and do **not** import
correctly; `predict_HE_unet.py` now has the same flag, still defaulting to `um`
so its behaviour is unchanged.

---

## Reproducibility

`set_seed` seeds python/numpy/torch(+CUDA) and the DataLoader uses a seeded
`Generator`, so **CPU runs are bit-exact at a fixed `--seed` with
`--num-workers 0`** (the default; the dataset caches whole-donor arrays, which
every forked worker would duplicate). GPU is reproducible only up to cuDNN/AMP.
`cv_summary.json` embeds the full per-epoch history (loss / IoU / lr) so under-
vs over-training is inspectable after the fact.

---

## What is not in this pipeline

* **SmallVessel segmentation.** Labelled at native resolution; sub-pixel at
  2 µm/px. Needs its own model at its own scale.
* **Cross-platform harmonisation.** The two 0.374 µm/px CODEX-panel images
  (`1901HBMP004`, `HDL073_PC29`) are exported onto the same 2 µm/px grid as the
  Phenocycler slides, but they are smaller sections of different image quality.
  H7's platform-sensitivity split is the precedent if this turns out to matter.
* **A validation notebook.** An `H20_unet_validation.ipynb` comparing manual,
  H&E-U-Net and DAPI-U-Net follicle metrics on the four donors present in both
  cohorts (HDL055/063/073/086), stratified by rs3184504, is the natural next
  step. `report_DAPI_concordance.py` covers the single-cohort half of it.
* **Groovy in CI.** No QuPath on the runner. The parity contract is covered by
  the transcription tests and by `check_dapi_parity.py` on real exports.
