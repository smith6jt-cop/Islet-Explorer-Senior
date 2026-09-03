# Spatial-method benchmark (BANKSY · UTAG · KODAMA · MISTy)

Compares the **relative performance** of four spatial-transcriptomics methods —
plus the in-pipeline niche baseline — on this project's Xenium 5k-panel samples
(`0041323` pancreas, `0041435` thymus). All four methods are run on **each**
sample; the primary comparison is **between methods, within each tissue**.

"Performance" = **computational** (runtime / peak RSS / scalability) **and
analytical quality** (domain coherence, cross-method agreement, robustness,
biological recovery; MISTy on its own R²/importance terms).

## Two tracks + bridge

| Method | Repo | Output | Track |
|---|---|---|---|
| BANKSY | prabhakarlab/Banksy_py (`pybanksy`) | spatial domains (λ=0.8) | domain |
| UTAG | ElementoLab/utag | spatial domains | domain |
| KODAMA | tkcaccia/KODAMA v3 | embedding → KMeans domains | domain |
| baseline | stage-03 niche (recomputed) | composition-KMeans domains | domain |
| MISTy | saezlab/mistyR | multiview predictive model | interaction |

- **Domain track** (BANKSY/UTAG/KODAMA/baseline): head-to-head, granularity-matched to ~12 domains.
- **Interaction track** (MISTy): per-target `intra/multi/gain.R2`, view contributions, predictor→target importances — never forced into a domain ranking.
- **Bridge**: Spearman(per-marker `gain.R2`, per-marker Moran's I) links the two tracks.

**Reporting stance (Tier-C):** objective metrics are tabulated; the **user**
ranks. No proportion-plausibility verdict, no "winner" cell. Every domain-quality
metric carries a bias note (see the generated report).

## Environments (production env untouched)

`setup_envs.sh` builds two isolated environments — it never modifies
`xenium_analysis` (whose scanpy 1.12 / rapids pins are load-bearing):

- **`bench_py`** conda env — BANKSY + UTAG (UTAG hard-pins `scanpy<1.12`).
  Their matplotlib needs a newer `libstdc++`; `env_bench_py.sh` sets the
  required `LD_PRELOAD`. `setuptools<81` is pinned (squidpy→xarray_schema needs
  `pkg_resources`).
- **R packages in `/blue/maigan/smith6jt/R_libs`** — `mistyR` (Bioconductor) +
  `KODAMA` (GitHub), installed against the existing R 4.5.3. An R user-library
  cannot touch Python site-packages, so the pins stay safe.

```bash
bash scripts/benchmark/setup_envs.sh          # both (idempotent)
bash scripts/benchmark/setup_envs.sh --py-only # just bench_py
```

## Run

```bash
# headline: 50k contiguous crop, both samples, all methods + eval
STEPS="prep methods eval" SAMPLES="0041323 0041435" VARIANTS="crop50k" \
  bash scripts/benchmark/run_benchmark.sh

# one method / one sample (e.g. after an edit)
SAMPLES=0041435 VARIANTS=crop50k STEPS=utag bash scripts/benchmark/run_benchmark.sh

# scalability sweep (then eval reads the curve)
VARIANTS="crop25k crop50k crop100k crop200k" STEPS="prep methods" bash scripts/benchmark/run_benchmark.sh
SWEEP_VARIANTS="crop25k crop50k crop100k crop200k" STEPS=eval VARIANTS=crop50k bash scripts/benchmark/run_benchmark.sh
```

SLURM: `sbatch scripts/benchmark/slurm/bench.sh` (all-in-one) or the
prep→array→eval chain documented in `slurm/bench_array.sh`.

### Variants
`crop<N>k` = compact contiguous crop (N nearest the spatial centre — preserves
the micron neighbourhood graph); `thin<N>k` = random N (sampling control);
`full` = all cells. Headline = `crop50k`. `BENCH_TIMEOUT=<s>` caps a run and
records `status="capped"` (no silent drops).

## Pipeline

```
prep_benchmark_input.py  -> data/processed/benchmark/{s}/{s}_bench_{v}_{input.h5ad,expr.mtx,coords.csv,meta.csv,barcodes.csv,genes.csv,misty_markers.txt,cropinfo.json}
methods/run_*.{py,R}     -> results/{method}_{v}_labels.csv  (+ misty_{v}/ + tidy CSVs)
profile_wrap.sh          -> results/{method}_{v}_metrics.json  (wall_s, peak_rss_mb, n_domains, status)
evaluate.py              -> {s}_benchmark_scorecard_{v}.csv, {s}_benchmark_report_{v}.md, figures/{s}/benchmark/*.png
```

All methods consume ONE identical input per `(sample, variant)` (genes =
`var['embedding_panel']`, the curated lineage-discriminating panel; configurable
via `GENE_PANEL`). MISTy uses the tissue `MARKER_DICT` union subset.

## Outputs
- **Scorecard** `data/processed/benchmark/{s}/{s}_benchmark_scorecard_{v}.csv`
- **Report** `…/{s}_benchmark_report_{v}.md` (scorecard + bias notes + agreement + MISTy + scalability)
- **Figures** `figures/{s}/benchmark/` (spatial maps, quality bars, AMI/ARI heatmaps, enrichment, MISTy plots, scalability, stability)

## Notes / limitations
- **KODAMA is the slowest method** by far. It runs with reduced, documented
  settings (`M=50`, `landmarks=2000`) for tractability; full `M=100`/
  `landmarks=10000` is hours on ≥25k cells. Override with `--M` / `--landmarks`.
- All four methods are **CPU**; VRAM is N/A (the only optional GPU use is Moran's I).
- **No ground-truth domains** → metrics are reference-free or soft-reference
  (vs pipeline celltypes); no accuracy claims.
- **Stability** (seed ARI) populates only if seeded runs exist under
  `results/stability/` (optional extra runs).
- Cross-tissue is descriptive only (n=2 tissues; no significance test).
