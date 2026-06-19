# HANDOVER — RESTORE per-image intensity normalization

> Ephemeral cross-instance handover (written 2026-06-18). **Delete after reading; do not commit.**
> Status: **COMPLETE & validated.** Everything below is on disk and working.

## TL;DR

Built a new standalone analysis that applies **RESTORE** (per-image autofluorescence
auto-thresholding for multiplexed imaging; Chang Lab/OHSU) to the 59-plex protein
single-cell data (`data/cells/`, 23.3M cells, 15 donor images) for **9 mutually-exclusive
markers** spanning the broad pancreas cell types. Produces per-image thresholds, per-cell
gating booleans, and a normalized intensity layer. The result independently recovers the
two T1D hallmarks (β-cell loss, T-cell infiltration), which validates it.

The work is done. This doc exists so you can continue/extend without re-deriving anything.

## What the user asked for (chronological)

1. "Use https://github.com/smith6jt-cop/RESTORE with the mutually-exclusive markers CD3e,
   CD99, CD163, CD20." → built + ran (chose SSC model, see below).
2. "CD99 needs to be dropped and repeated without it. Suggest mutually-exclusive pairs for
   the broadest pancreas cell types." → dropped CD99, data-mined marker choices, added 6
   broad-lineage pairs, re-ran. **This is the current/final state.**

## The deliverable

**Main script:** `scripts/senior/restore_normalize.py` (untracked/new). Run in conda env
`scvi-env`. Self-contained: builds the RESTORE input, runs it, flattens thresholds, applies
to all cells, writes QC. Key CLI flags: `--model {SSC,GMM,KMeans}` (default **SSC**),
`--robust/--no-robust` (default on), `--robust-factor` (3.0), `--subsample` (15000),
`--reuse-threshs`, `--limit-scenes N`, `--skip-apply`, `--marker-pairs 'A:B,C:D'`.

**Final marker pairs** `[target ← mutually-exclusive reference]` (in `DEFAULT_MARKER_PAIRS`):
```
CD3e            ← CD163              # T cells       (immune refs kept immune-vs-immune)
CD20            ← CD3e              # B cells
CD163           ← CD3e              # macrophages
Pan_Cytokeratin ← Vimentin         # exocrine epithelium  (acinar/ductal)
Vimentin        ← Pan_Cytokeratin  # mesenchyme/stroma
INS             ← Pan_Cytokeratin  # endocrine beta
GCG             ← Pan_Cytokeratin  # endocrine alpha
CD31            ← Pan_Cytokeratin  # endothelium
SMA             ← Pan_Cytokeratin  # smooth muscle / pericyte
```

## Environment (CRITICAL — RESTORE doesn't install cleanly)

- Env: `conda activate scvi-env` (python 3.13; has numpy/pandas/pyarrow/scipy/sklearn 1.7/
  matplotlib/anndata).
- **RESTORE is vendored, NOT pip-installed** (its `setup.py` is a broken placeholder):
  ```bash
  git clone https://github.com/smith6jt-cop/RESTORE.git external/RESTORE   # pinned @ 38df59b
  pip install spams-bin     # SSC model needs `spams`; spams-bin 2.6.13 cp313 wheel installed OK
  ```
  `external/` is gitignored (added this session). Already cloned & on disk.
- RESTORE's `RESTORE.py` imports `holoviews`, `selenium`, `tqdm.notebook` at module load AND
  builds `hv.Curve/Scatter/Overlay` objects *during* threshold computation (not just when
  saving figs). The script injects lightweight **stub modules** (see `_install_headless_stubs()`)
  so the vendored source runs headless **unmodified** with `save_figs=False`. It also prepends
  `external/RESTORE/python_code` to `sys.path` (resolves both `RESTORE` and bundled `ssc`).
  → Do NOT edit the vendored source; the stubs handle everything.

## How RESTORE works (so you don't have to read its source)

For target marker T paired with reference R: R⁺ cells are (by biology) T-negative, so their
T-intensity is background. Per image, fit a 2-cluster model on raw 2D `[T, R]` (after an
`idx_select` prefilter, abs floor `>50`), pick the R-dominant (T-negative) cluster, set
`threshold = mean + 3·σ` of T in that cluster. Three models computed (KMeans/GMM/SSC),
all saved. Output dict: `threshs[batch][scene][target][model]`. We map `scene=donor_id`,
`batch=0` (single batch; `threshs[0][...] == threshs['global'][...]`). **Operates on RAW
MFI — do NOT log-transform input.**

## Key decisions + rationale (don't relitigate)

- **SSC model chosen** (over RESTORE's usual GMM). On this panel GMM/KMeans place abundant
  markers' thresholds at/above the 99.9th intensity pct → ~0% positive (implausible). SSC
  lands at the background→signal shoulder. Verified by comparing thresholds to raw
  percentiles and per-image positive fractions. All 3 models still saved in the CSV.
- **Robustness guard (default on).** A few per-image SSC thresholds are degenerate (abundant
  Pan_CK/Vimentin overshoot to 3000–6400 vs ~850 median; INS undershoots). `--robust`
  imputes any threshold `>3×` or `<1/3×` the cohort median with the cohort median, logged
  per imputation. Raw per-model thresholds in `restore_thresholds.csv` are left untouched.
  24/135 imputed in the 9-marker run.
- **Marker selection is data-driven** (sampled 150k cells; checked expression + Spearman
  exclusivity). Pan_Cytokeratin = universal reference (~91% of cells = acinar/ductal,
  cleanly negative elsewhere). Epithelial↔mesenchymal (Pan_CK↔Vimentin) corr ≈ −0.5 = the
  broadest exclusivity. **Rejected** as poor exclusive anchors: CD99 (96%+ detectable, not
  lineage-specific — this is why it was dropped), CD56, CD34, EpCAM, SST, IAPP, CD11b, B3TUBB.
- **Immune references kept immune-vs-immune** (user choice), not switched to Pan_CK.

## Outputs (all under gitignored `data/`, already written)

- `data/restore_thresholds.csv` — tidy, all 3 models × 15 images × 9 markers = **405 rows**
  (`image,marker,model,threshold,chosen`).
- `data/restore_gated/donor_id=*/data_0.parquet` — Hive-partitioned, **29 cols**: `object_id,
  donor_id` + for each of 9 markers `{<m>_pos (bool), <m>_norm (=raw/threshold, threshold→1.0),
  <m>_log2r}`. ~2.6 GB, 15 partitions, 23.3M rows total.
- `data/restore/threshs.pkl` (native RESTORE), `data/restore/positive_fractions.csv`
  (has an `imputed` flag column), `data/restore/qc/*.png` (17: 15 per-image histograms +
  `threshold_heatmap.png` + `positive_fractions.png`).

## How to run / re-run

```bash
conda activate scvi-env
cd /home/smith6jt/IO60panc2nd/Islet-Explorer-Senior

python scripts/senior/restore_normalize.py                                # full 9-marker run (~5 min)
python scripts/senior/restore_normalize.py --reuse-threshs --no-robust    # re-apply, pure thresholds (fast, ~30s)
python scripts/senior/restore_normalize.py --limit-scenes 1 --skip-apply  # single-image dry run
python scripts/senior/restore_normalize.py --marker-pairs 'CD3e:CD163,INS:Pan_Cytokeratin'  # custom
```
Thresholds are deterministic (seed=0) → `--reuse-threshs` is safe for re-applying with
different `--model`/`--robust` options without recomputing.

## Validation (the result that matters)

Per-image RESTORE fractions vs donor `Donor Status` (from `singlecell_protein.h5ad.obs`),
recovered WITHOUT seeing disease labels:

| Status | INS⁺ (β cells) | CD3e⁺ (T cells) |
|--------|----------------|-----------------|
| ND     | 3.39%          | 0.82%           |
| Aab+   | 2.30%          | 1.10%           |
| T1D    | 0.48%          | 3.02%           |

β-cell loss + T-cell infiltration across ND→T1D = exactly what this study tracks. T1D
images 6380/6579 show 0.0% INS (genuine beta loss, not artifact).

## Caveats / known limitations

- **Pan_Cytokeratin⁺ ≈ 20–35%** (conservative for the exocrine majority). `mean + 3σ`
  flags only high-confidence epithelial cells. RESTORE hardcodes `sigma_weight=3`; if you
  want a fuller acinar compartment, expose/lower it (it's in `external/RESTORE/.../RESTORE.py`
  `get_GMM_thresh`/`get_clustering_thresh`; would need a small patch or a `--sigma-weight`
  passthrough, currently NOT plumbed through).
- INS thresholds vary biologically (beta loss) → robust imputation slightly raises them in
  T1D images but does not erase the signal. Use `--no-robust` for pure per-image values.
- No existing cell-type labels to cross-check against — **Senior phenotyping is still staged**
  (`singlecell_protein.h5ad` has no phenotype obs column). Natural next step: feed these
  `*_pos`/`*_norm` columns into that phenotyping.

## Repo state / uncommitted changes (NOT committed — user hasn't asked)

```
 M .gitignore                          # added `external/` (vendored, gitignored)
 M scripts/CLAUDE.md                   # added "RESTORE per-image intensity normalization" section
 M scripts/senior/cluster_annotate.py  # PRE-EXISTING (was modified before this session; NOT mine — leave it)
?? scripts/senior/restore_normalize.py # the new script
```
Plus untracked gitignored artifacts: `external/RESTORE/`, `data/restore/`,
`data/restore_thresholds.csv`, `data/restore_gated/`.

If asked to commit: branch off `main`, commit only `scripts/senior/restore_normalize.py`,
`scripts/CLAUDE.md`, `.gitignore`. **Do NOT** touch `cluster_annotate.py` (unrelated WIP)
or commit anything under `data/`/`external/`.

## Possible next steps the user may want

1. Commit the script + doc/gitignore edits (see above).
2. Plumb a `--sigma-weight` flag through to RESTORE (needs a vendored patch) to make
   Pan_CK less conservative.
3. Feed `data/restore_gated/` into Senior phenotyping (its real downstream purpose).
4. Aggregate per-cell gating to per-islet composition (joins on `object_id` / islet ids in
   `data/cells/`).

Original plan file (this instance): `/home/smith6jt/.claude/plans/in-a-new-analysis-harmonic-axolotl.md`
Full per-script reference: `scripts/CLAUDE.md` § "RESTORE per-image intensity normalization (Jun 2026)".
