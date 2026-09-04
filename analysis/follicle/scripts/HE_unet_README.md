# H&E Follicle U-Net Pipeline

Automated follicle segmentation on H&E whole-slide images of human spleen.
Trains on the same per-follicle polygons that the H&E QuPath project already
exports, and emits output in the schema `H18_HE_follicles.ipynb` already reads.

> **Single-channel hematoxylin (H) input.** The model trains and infers on the
> deconvolved hematoxylin concentration channel (`in_channels=1`), **not** RGB.
> Follicle detection is a nuclear-density problem, so eosin and inter-donor
> colour are nuisance variation. The upstream Macenko RGB normalization is kept;
> the H channel is derived from the normalized RGB (high H → dark). The QuPath
> export writes the H image (plus an RGB-derived tissue mask) directly; the
> Python `rgb_to_h` is the inference-side counterpart and **must stay
> numerically identical to the export** (see *Stain deconvolution & parity*).

## Why this exists

The fluorescent cohort (`AnnotationsFinal.csv`) gets follicles from the
`SpleenRegions2` pixel classifier (`scripts/Workflow.groovy`). The H&E cohort
(`27_donor_spleen_measurements1.csv`) does not — every follicle in that file
was drawn by hand. This pipeline replaces the manual step with a U-Net that
matches the existing schema, so it slots into `H18` and the H&E cohort of
`H19_meta_analysis.ipynb` without changing the downstream notebooks.

## Files

| File | Purpose |
|---|---|
| `export_HE_training_pairs.groovy` | QuPath script: per H&E entry, deconvolves the normalized RGB to a single-channel hematoxylin PNG + an RGB-derived tissue mask PNG (`data/tissue/`) + matching binary Follicle mask PNG + JSON sidecar (with the `h_transform` block) into `data/`. Runs once. |
| `HE_unet_utils.py` | Library: `rgb_to_h` (inference-side RGB→H, parity-locked to the export), `FolliclePatchDataset`, `build_unet` (`in_channels=1`, `pretrained=` toggle), `bce_dice_loss`, `predict_wsi` (single-channel H in; tiled inference with cosine blending + 4-way TTA), `extract_follicle_records` (mask → instance polygons via watershed), `instance_metrics` (count/area vs GT), `write_geojson`, `donor_kfold` (genotype-stratified), `set_seed`/`seed_worker`, `check_mask_integrity`, `discover_samples(data_dir=…)`. |
| `train_HE_unet.py` | Donor-level k-fold CV training. Saves `analysis/HE_unet_models/fold_{i}.pt` + `cv_summary.json` (with full per-epoch history + best-checkpoint count/area vs GT) + `pooled_metrics.json`. Supports `--smoke`, `--patience` early stopping, `--lr-tmax`, `--data-dir`/`--out-dir`, `--no-pretrained`. |
| `predict_HE_unet.py` | WSI inference across one or more fold checkpoints. Reads the (normalized) slide as RGB, converts to H once at the boundary via `rgb_to_h`, then runs the model. Writes per-image GeoJSON to `analysis/HE_unet_predictions/` and a combined `Measurements/HE_follicles_unet.csv`. |
| `make_synthetic_he_data.py` | Generates synthetic (single-channel H image, mask, tissue, meta) tuples — the single source of truth for the pytest suite and the CI end-to-end smoke-train step. |
| `requirements_HE_unet.txt` | Pip dependencies for the whole pipeline (note: `scikit-image>=0.26`). |
| `../tests/test_HE_unet_utils.py` | CPU-only pytest suite. Runs in CI, alongside an end-to-end `--smoke` train on synthetic data. |

## Quickstart

### 1. Install

```bash
python -m venv .venv && source .venv/bin/activate
pip install -r scripts/requirements_HE_unet.txt
```

For GPU training, replace the `torch` line in `requirements_HE_unet.txt` with
the CUDA wheel matching your driver:

```bash
pip install torch --index-url https://download.pytorch.org/whl/cu121
```

`openslide-python` is a thin wrapper over the OpenSlide C library — on Ubuntu
you also need `apt install libopenslide0`. macOS: `brew install openslide`.
The pipeline falls back to `tifffile` for plain TIFF/OME-TIFF inputs when
OpenSlide is unavailable.

### 2. Export training pairs from QuPath

Open the H&E QuPath project that contains the manual follicle annotations,
then run from a terminal using QuPath's bundled command-line launcher
(`bin/QuPath` inside the install on Linux/macOS, `QuPath.exe` on Windows):

```bash
QUPATH="$HOME/QuPath-v0.7.0-Linux/QuPath/bin/QuPath"
"$QUPATH" script --project path/to/HE_project.qpproj scripts/export_HE_training_pairs.groovy
```

(No `--save` is needed here — the export only reads annotations and writes
the PNG/JSON pairs to disk.)

Result: `data/{images,masks,tissue,meta}/<donor>.{png,png,png,json}`
for every project entry that has at least one `Follicle` annotation, where
`images/<donor>.png` is the **single-channel hematoxylin** image (high H =
dark) and `tissue/<donor>.png` is the RGB-derived tissue mask
(`mean(R,G,B) < 220`, so red pulp is correctly tissue). All donors land on a
common 2 µm/px grid regardless of native scanner resolution. The export must
run on the **same normalized images** that inference will see.

The directory is already in `.gitignore` so the (sometimes hundreds of MB
of) PNGs stay out of the repo. **Existing RGB exports are the wrong format —
regenerate `data/` after switching to the H pipeline (see *Migration*).**

### 3. Train

```bash
# Fast smoke test — synthetic data, CPU, no pretrained weights, ~seconds.
# A crash-check of the whole pipeline, NOT a quality run.
python scripts/make_synthetic_he_data.py --out-dir /tmp/he_synth --n-donors 3
python scripts/train_HE_unet.py --smoke --data-dir /tmp/he_synth

# Single-fold quality estimate on the real data — train to a plateau.
python scripts/train_HE_unet.py --folds 1 --epochs 60 --patience 12

# Full genotype-stratified 5-fold CV — the reportable run.
python scripts/train_HE_unet.py --folds 5 --epochs 50 --batch-size 16
```

A single-fold run holds out just one donor, so its IoU is a sanity check, not
the headline number. The success bar — per-donor follicle count + mean-area
concordance with the manual annotations (Spearman ρ > 0.9) — needs the full
5-fold run; pooled count/area Spearman is written to `pooled_metrics.json`.

Outputs: `analysis/HE_unet_models/fold_{0..4}.pt` plus `cv_summary.json` listing,
per fold, the best validation IoU, the **full per-epoch history**
(loss / val IoU / lr), and best-checkpoint follicle count + mean area vs.
ground truth.

Useful flags:

- `--smoke` — fast crash-check preset: tiny epochs/patches/batch, CPU, no
  pretrained weights (no network). Pair with `--data-dir`/`--out-dir`.
- `--data-dir <dir>` / `--out-dir <dir>` — read training pairs from / write
  checkpoints to a custom location (defaults: repo `data/` and
  `analysis/HE_unet_models/`). Keeps smoke/CI runs off the git-tracked dirs.
- `--patience N` — stop a fold after N epochs with no val-IoU improvement
  (0 = disabled, the default; `--smoke` uses 1).
- `--lr-tmax N` — cosine-annealing period in epochs (default: `--epochs`). Set
  it to your expected effective budget so the LR doesn't bottom out early when
  you also pass `--patience`.
- `--no-pretrained` — random-init encoder, skipping the ImageNet download.
- `--encoder efficientnet-b3` (default) — any encoder name from
  `segmentation_models_pytorch` works. Falls back to a plain-PyTorch U-Net
  if SMP isn't installed.
- `--genotype-xlsx <path>` — used to stratify the donor-level folds so the
  9/9/9 rs3184504 balance survives each split. Defaults to
  `Spleen_rs3184504_Genotypes.xlsx` at the repo root.
- `--num-workers 0` (default) — the dataset caches whole-donor arrays in
  RAM, so each forked DataLoader worker would duplicate that cache. Bump
  only if you have confirmed memory headroom.

### 4. Predict over whole slides

```bash
python scripts/predict_HE_unet.py \
    --wsi-dir /path/to/HE_wsis \
    --models analysis/HE_unet_models/fold_*.pt
```

Outputs:

- `analysis/HE_unet_predictions/<wsi>.geojson` — drop into QuPath via
  *Objects → Import objects* to QC visually.
- `analysis/HE_unet_predictions/<wsi>.meta.json` — pixel size, follicle
  count, threshold, ensemble size.
- `Measurements/HE_follicles_unet.csv` — combined records using the same
  schema as `27_donor_spleen_measurements1.csv`, ready for H18.

Useful flags:

- `--threshold 0.5` — probability cutoff (default).
- `--min-area-um2 10000` — drops connected components below this (default).
- `--no-tta` — disables 4-way test-time augmentation. ~4× faster, ~1–2 IoU worse.
- `--no-watershed` — disables touching-follicle splitting. Use if the
  watershed seeds are over-segmenting on your slides.

### 5. Drop into H18

Change one line in `analysis/H18_HE_follicles.ipynb`:

```python
he = pd.read_csv(du.DATA_ROOT / "Measurements" / "HE_follicles_unet.csv")
```

If U-Net-vs-manual concordance is high (Spearman ρ > 0.9 on per-donor count
and mean area), H18 conclusions reproduce. If not, the discrepancy is the
finding worth investigating.

## Design choices and tradeoffs

### Working resolution: 2 µm/px

Follicles are 100–1000 µm across, so a 512² patch covers ~1 mm. At native
0.5 µm/px we'd waste capacity on hematoxylin texture; at 4 µm/px we'd lose
follicle boundary precision. 2 µm/px is the sweet spot for this structure.

### Single-channel hematoxylin input & export↔Python parity

The model input is the deconvolved hematoxylin (H) concentration channel, a
single channel (`in_channels=1`). `smp.Unet` adapts the ImageNet-pretrained
first convolution to one input channel by summing the RGB filter weights, so the
encoder stays pretrained. Eosin and inter-donor colour are nuisance variation
for a nuclear-density problem, and the upstream Macenko RGB normalization
already handles cross-donor stain intensity — a fixed reference H&E basis is
sufficient downstream.

The RGB→H transform exists in two places that **must stay numerically
identical**:

- `export_HE_training_pairs.groovy` — produces the on-disk training H images.
- `rgb_to_h` in `HE_unet_utils.py` — the inference-side counterpart, applied to
  whole-slide RGB in `predict_HE_unet.py` immediately after the slide is read.

Both use the reference H&E OD unit vectors
`H = [0.5626, 0.7201, 0.4062]`, `E = [0.2159, 0.8012, 0.5581]` (normalized),
`OD = -log10(max(channel, 1) / 255)`, `cH = w·OD` where `w` is the first row of
`(MMᵀ)⁻¹M` for `M = [Ĥ; Ê]`, `cH` clipped to `[0, 2.0]`, and intensity
`= floor(255·(1 − cH/2.0))` (high H → dark). **Any change to these constants
must be made in both files in the same commit.** Verify byte-for-byte parity:

```python
import numpy as np, sys; sys.path.insert(0, "scripts")
from HE_unet_utils import rgb_to_h
H = np.array([0.5626,0.7201,0.4062]); E = np.array([0.2159,0.8012,0.5581])
Hn, En = H/np.linalg.norm(H), E/np.linalg.norm(E)
hh, he, ee = Hn@Hn, Hn@En, En@En
det = hh*ee - he*he; a00, a01 = ee/det, -he/det
w = a00*Hn + a01*En
lut = np.array([-np.log10(max(v,1)/255.0) for v in range(256)])
def export_h(rgb):
    flat = rgb.reshape(-1,3); out = np.empty(len(flat), np.uint8)
    for i,(r,g,b) in enumerate(flat):
        cH = min(max(w[0]*lut[r]+w[1]*lut[g]+w[2]*lut[b], 0.0), 2.0)
        out[i] = int(min(max(255.0*(1-cH/2.0), 0.0), 255.0))
    return out.reshape(rgb.shape[:2])
rgb = np.random.default_rng(0).integers(0,256,(64,64,3)).astype(np.uint8)
assert np.array_equal(rgb_to_h(rgb), export_h(rgb)), "PARITY BROKEN"
print("export/python H parity OK")
```

CI does not run Groovy, so this parity check is **manual** — run it whenever the
export or `rgb_to_h` changes.

### Tissue mask

The export writes `data/tissue/<donor>.png` from the RGB as `mean(R,G,B) < 220`,
so red pulp (low H, near-white in the H channel) is still counted as tissue when
the patch sampler picks true-negative tiles. The dataset reads this mask; it
only falls back to deriving tissue from the H image alone (high H = dark) when
the PNG is absent — which under-counts red pulp, hence the dedicated mask.

### Donor-level splits, not patch-level

Two patches from the same WSI are not independent samples. `donor_kfold`
splits by donor, stratified on rs3184504 so each fold sees the 9/9/9
genotype balance. Patch-level CV would inflate apparent IoU by ~10 points
and tell you nothing about generalization to new patients.

### EfficientNet-B3 encoder, BCE + Dice loss

Pretrained ImageNet encoders consistently outperform from-scratch at this
data scale (n=27 donors). EfficientNet-B3 beats ResNet-34 on most histology
benchmarks for similar VRAM. Loss is `0.5 * BCE + 0.5 * Dice` because pure
BCE is unstable on ~5–10% positive pixels and pure Dice is noisy at small
batch sizes.

### Touching-follicle separation via distance-transform watershed

Biological follicles can touch. After thresholding the probability map we
run watershed seeded on local maxima of the distance transform, with a
~50 µm minimum centre-to-centre separation. When `peak_local_max` returns
no seeds (single tiny blob, flat distance plateau) the code falls back to
`skimage.measure.label` so no instances are silently dropped.

### Output schema

`Measurements/HE_follicles_unet.csv` columns:

| Column | Meaning |
|---|---|
| `Image` | WSI filename |
| `Object ID` | `unet-<image>-<label>` |
| `Object type` | `Annotation` |
| `Name` | empty |
| `Classification` | `Follicle` |
| `Parent` | `Root object (Image)` |
| `ROI` | `Polygon` |
| `Centroid X µm`, `Centroid Y µm` | in µm |
| `Area µm^2` | in µm² |
| `Perimeter µm` | in µm |

Identical to `27_donor_spleen_measurements1.csv` so H18 needs no schema
changes.

## Migration (switching an existing RGB setup to H)

These are operational steps, not code edits, and are required for correctness:

1. **Regenerate `data/`.** Any existing `data/images/*.png` are RGB and are now
   the wrong format/contract. Re-run `export_HE_training_pairs.groovy` over the
   QuPath project to produce H images and tissue masks. (`data/` is gitignored,
   so this is a local action.)
2. **Existing checkpoints are invalid.** A model trained on RGB has a 3-channel
   first conv and will not load into the `in_channels=1` model. Retrain from
   scratch (the encoder first conv is re-derived for one input channel anyway).
3. **Point inference at normalized slides.** `predict_HE_unet.py --wsi-dir` must
   point at the **same normalized** whole slides the training images came from,
   never raw scans — the model's H statistics assume normalized input.

## What's not in this pipeline

- **Stain normalization.** The upstream Macenko RGB normalization is kept and
  lives outside this repo; source filenames already say `_normalized.ome.tif`.
  The export derives the H channel from that normalized RGB, so cross-donor
  stain intensity is handled before deconvolution. The related (out-of-repo)
  stain-region scorer was switched to an H-bimodality criterion consistent with
  the H-only decision; no action here unless that script is added to the repo.
- **PALS (T-zone) segmentation.** Single-class model only. Extending to a
  2-channel head is straightforward but needs PALS annotations on the H&E
  cohort first.
- **U-Net-vs-manual concordance notebook.** `H20_HE_unet_validation.ipynb`
  is a planned follow-up once a trained model exists. It should report IoU
  and per-donor concordance separately by rs3184504 — a model that performs
  worse on T/T donors would silently bias H18.
- **An H&E QuPath project.** This pipeline assumes one already exists with
  manual Follicle annotations in place. The 27-donor CSV in the repo was
  exported from one; the project itself lives outside this repo.

## CI

`.github/workflows/he_unet_ci.yml` runs `tests/test_HE_unet_utils.py` on
Python 3.11 and 3.12 with CPU-only torch on every push or PR that touches
the pipeline. The suite covers post-processing geometry, watershed edge
cases (the ones Copilot's review flagged), donor-level fold partitioning,
loss numerics, reproducibility of the seeded data pipeline, the mask-integrity
preflight, the deployment-aligned instance metrics, and a tiny forward pass
through the fallback U-Net. A final step then runs `train_HE_unet.py --smoke`
end-to-end on freshly generated synthetic data, so the training entrypoint
itself — not just the unit-tested helpers — is exercised on every change.

## Reproducibility

- `--seed` (default 0) seeds python/numpy/torch (+CUDA) via `set_seed`, and the
  training DataLoader gets a seeded `torch.Generator` — so both the data order
  and the randomly-initialised decoder are fixed run-to-run.
- **CPU runs (`--device cpu`) are bit-exact across reruns at the same seed**
  (with the default `--num-workers 0`).
- GPU runs are reproducible only up to cuDNN/AMP nondeterminism — two reruns may
  differ in the 3rd–4th decimal. `set_seed` sets `cudnn.deterministic=True`; we
  deliberately do **not** call `torch.use_deterministic_algorithms`, which needs
  `CUBLAS_WORKSPACE_CONFIG` set before the process starts and can hard-error on
  some encoder upsampling kernels.
- Each fold records its full per-epoch history in `cv_summary.json`, so the
  loss / val-IoU curve — and whether a run was under- or over-trained — is
  inspectable afterward, not just the single best value.
- Checkpoints record the encoder name, fold index, and best per-donor IoU
  so prediction is fully traceable to the training run that produced it.
