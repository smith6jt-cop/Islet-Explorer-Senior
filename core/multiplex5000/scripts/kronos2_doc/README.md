# `scripts/kronos2_doc/` — KRONOS2 run the way CORAL documents it

Maintainer instruction, 2026-09-07: *"run KRONOS2 according to the documentation and let me verify the
results."* This lane is **not** a bake-off arm. There are no ablations, no null gate, no verdict and no
release: it is CORAL's own documented pipeline, run as documented, so its clusters and their marker
signatures can be read directly. The arm `broadtype/arms/kronos2_coral.py` and its records are untouched.

```bash
cd core/multiplex5000
../../.venv/bin/python scripts/kronos2_doc/run.py --donors 6436 6442
../../.venv/bin/python scripts/report/run_report.py --section kronos2doc     # K01-K03
```

## What is run, and where it comes from

| step | source | what runs |
|---|---|---|
| ingest | Tutorial 1 | `convert_to_canonical` per frame window |
| tissue | Tutorial 2 (a stated prerequisite of Tutorial 3) | `detect_tissue(OtsuTissueSegmenter())`, the CLI's default method |
| cells | Tutorial 3 step 1, Approach B | `import_cell_mask` — ids preserved, mask authoritative, not tissue-restricted |
| patches | Tutorial 3 step 2 | `PatchConfig(patch_size=64, mode="cell_centered")` = `coral patch --mode cell --patch-size 64` |
| features | Tutorial 3 step 3 | `mean_marker` baseline **and** `KRONOS2`, both on the panel subset |
| clusters | Tutorial 4 | `sc.pp.neighbors(n_neighbors=15, use_rep="X")` → `sc.tl.leiden(resolution=0.5, flavor="igraph", n_iterations=2, directed=False)` → `sc.tl.umap`, all at seed 0 |

Three environments, each owning its step: `.venv` cuts the windows out of the qptiff,
`external/.venv-kronos` runs CORAL and KRONOS2, `scvi-env` runs the scanpy clustering (scanpy, igraph
and umap-learn are deliberately **not** CORAL dependencies — CORAL's scope ends at the feature matrix).

## Deviations, all recorded in `features_record.json`

- **Scope.** The tutorial ingests one slide. A PhenoCycler slide here is 23,760 × 24,960 px × 59 channels
  and `convert_to_canonical` reads it as a single array (~70 GB), so the unit is the frozen screen frame:
  sixteen 2048-px windows per donor, the same cell set every bake-off arm was judged on.
- **Batch size.** The tutorial gives two values — `--batch-size 1024` on its CLI line and `batch_size=16`
  in the Python on the same page. CORAL's source settles it: *"the gold-standard features were produced
  at batch 16 … reproduce the gold standard with `batch_size=16`"*, and batch invariance is
  GPU-dependent (drift up to ~1e-4). 16 is what runs. 1024 would not fit a 20 GB card at 29 markers × 64 px.
- **Patch size.** 64 px as documented. The tutorial notes this is "~24 µm here" at its slide's 0.37 µm/px;
  at 0.5078 µm/px the same 64 px is 32.5 µm. The documented number is kept, the micron size recorded.
- **Panel.** DAPI + the 28 declared identity markers, the analogue of the tutorial's "18 phenotypic
  markers — the panel that actually separates cell types — rather than all 49 channels". CORAL hard-errors
  on a selection without a nuclear stain, so DAPI is in it.
- **Novel markers.** Seven panel markers are outside KRONOS2's 268-marker vocabulary (Ker8_18, Keratin_5,
  Caveolin, TP63, GCG, INS, SST). CORAL's documented `prepare_cohort_stats` computes their (mean, std)
  from each donor's own tissue-masked windows and registers them; the compartment/family text rows reuse
  existing KRONOS2 categories and are recorded assumptions, not measurements.

## What this lane does NOT do

- **It assigns no cell types.** Tutorial 4 assigns none either — CORAL lists naming as future work
  ("the cluster × marker signature turns each unsupervised cluster into a phenotype") and states plainly
  that *"CORAL's job ends at the feature matrix"*. K02 is that signature; reading it is the maintainer's step.
- **Tutorial 5 is not available here.** It is a supervised linear probe on ground-truth `cell_labels.csv`.
  This project has no such labels it may train on — QuPath classifications are audit context only — and
  Tutorial 6 (LoRA) needs strictly more labels, not fewer.
- **Nothing here feeds a release.** No calls, no `evidence_positive`, no verdict; the registry latch is
  untouched.

## Worth knowing before reading the output

CORAL's own benchmark puts KRONOS2 barely ahead of plain mean-marker intensities for cell phenotyping —
balanced accuracy 0.720 vs 0.702 (cHL), 0.742 vs 0.737 and 0.796 vs 0.790 (DLBCL) — and says so directly:
*"a foundation model does not automatically win"*, with the gain concentrated in classes intensity alone
cannot resolve. That is why the mean-marker baseline is clustered and plotted beside KRONOS2 here.
Leiden resolution 0.5 is the tutorial's value with no stated justification and no stability criterion;
it is not a validated choice, and "higher → more, smaller clusters" is the only guidance given.
