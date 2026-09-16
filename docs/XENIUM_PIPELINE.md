# Xenium analysis pipeline — gDNA-corrected, per-donor

**Status:** Stages 1–2 implemented and validated. Stages 3–6 specified, not yet built.
**Scope:** Xenium Prime 5k, per donor. Nothing in this pipeline is cohort-wide.
**Applies to:** the Pheno-Xeno pair — 0080870 (donor 6533) and 0081366 (donor 112).
Both donors are **female** (confirmed by the maintainer).

---

## 0. Read this first: why the pipeline changed

These bundles are flagged for **high genomic DNA control signal**, and the flag changes what a
transcript count means. Every number below was read from the bundles, not assumed.

| `metrics_summary.csv` | 0041323 (what thresholds were tuned on) | 0080870 | 0081366 |
|---|---:|---:|---:|
| `adjusted_genomic_control_probe_rate` | 0.00454 | **0.261** | **0.306** |
| `genomic_control_probe_counts_per_control_per_cell` | 0.000325 | 0.02520 | 0.02571 |
| est. false-positive transcripts/cell (incl. genomic) | **1.66** | **126.04** | **128.56** |
| median transcripts/cell | 309 | 433 | 403 |

**~76× more gDNA-derived false positives per cell than the sample every constant was
calibrated against.** About 29% of a median cell's counts are background.

### The model, and why it is trustworthy

The vendor's model is that every gene probe picks up gDNA at the same rate a genomic control
probe does. It reproduces exactly:

```
vendor  genomic_control_probe_counts_per_control_per_cell = 0.0252021
        0.0252021 x 5001 genes                            = 126.04
        vendor states                                     = 126.04   <- 0.00% error
        mean genomic_control_counts/cell 0.5292 / 21 probes = 0.0252  <- same number
```

Panel composition (`gene_panel.json`): **5,001 genes, 21 `genomic_control`, 40
`negative_control`, 3,291 `deprecated_codeword`**. Genomic controls are named
`Intergenic_Region_*`.

**Independently confirmed** against genes that cannot be expressed in pancreas
(`analysis/xenium/panel_audit.csv`, 106 exclusions). These are real gene probes, so what they
pick up IS the gene-level background:

```
sex_chromosome          0.79%   0.32x prediction   <- Y genes, female donor: NO template
skin_cornification      2.56%   1.03x              <- clean negatives land on the prediction
nephron_specific        2.92%   1.17x
photoreceptor_retinal   3.51%   1.41x
spermatogenesis         3.62%   1.46x
hepatocyte              4.94%   1.99x   <- shared endoderm lineage: REAL low expression
oocyte_germline         5.04%   2.03x
embryonic_pluripotency 10.10%   4.06x   <- not background at all
```

Clean categories sit on the prediction; the gradient above them is biology. **Do not average
all 106 exclusions into a "background factor"** — doing so gives a spurious 1.8x.

### The ground-truth negative set

Both donors are female, so all six Y probes have **no template of any kind** — every count is
background. This is the validation target for any detection rule:

| | `TBL1Y` / `USP9Y` / `EIF1AY` | `DDX3Y` / `KDM5D` |
|---|---|---|
| raw count >= 1 | **0.21 – 0.35%** | 1.47 – 2.08% |

The floor group gives the non-gDNA background (optical + decoding) ~0.3%; the elevation in
`DDX3Y`/`KDM5D` is cross-hybridisation from their highly-expressed X paralogs `DDX3X`/`KDM5C`.
And `0.3% + 2.49% gDNA ~= 2.8%` reproduces the clean negatives.

**Any proposed detection predicate must put the Y probes at or below ~0.3%.**

### Raw counts vs lognorm — they are not interchangeable

Measured on the same zero-template Y probes:

| gate | `DDX3Y` false detection |
|---|---|
| raw count >= 1 | 1.88% |
| **lognorm >= log1p(1)** | **0.55%** |

`normalize_total` scales a single stray count below 1 in any cell of above-median depth, and
gDNA scales with cell size — so the gate bites hardest exactly where background is worst. The
repo's **raw-count floors are the vulnerable ones**; its lognorm family is comparatively safe.

---

## 1. The correction

```
lambda_c  = (genomic_control_counts_c + a * mu_c) / (n_probes + a)     a = n_probes = 21
mu_c      = mean(genomic_control_counts | area decile of c) / n_probes
corrected_counts_c = max( transcript_counts_c  -  lambda_c * n_genes , 0 )
```

`lambda_c` is the expected gDNA counts **per gene** in cell c. Shrinkage is needed because the
raw estimate is a Poisson draw with mean 0.53 over 21 probes and **61.8% of cells record zero
genomic counts** — those are not zero-background cells. Area is the right conditioning
variable: gDNA per 100 µm² is roughly flat (0.98 → 0.77 across quintiles) while gDNA per cell
rises 4x (0.212 → 0.805).

**Three limits, all load-bearing:**
1. It subtracts an **expectation, not identified reads**. Valid for thresholds and aggregates;
   it cannot clean an individual cell's expression vector.
2. It inherits the vendor's uniformity assumption. Genes in high-copy or repeat-adjacent loci
   will exceed it.
3. `max(..., 0)` clips the 2.8% / 3.5% of cells whose expected background exceeds their entire
   transcript count. Flag them; do not silently zero them.

For **per-gene** decisions use `lambda_c` directly as a Poisson rate — `P(K >= k | lambda_c)` —
not multiplied by 5,001.

---

## 2. Stages

### Stage 1 — `scripts/xenium_celltype/gdna.py`  ✅ IMPLEMENTED, VALIDATED

```bash
cd core/phenocycler
python3 scripts/xenium_celltype/gdna.py \
    --bundle /path/to/output-XETG00298__0080870__panc__20260910__211917 \
    --out data/development/xenium_celltype/gdna_0080870.json
```

Reads `n_probes` from `gene_panel.json` (never hard-coded), computes `lambda_c`, and
**re-derives the vendor's headline false-positive number, failing loud above 5% disagreement.**
Verified: 126.04 vs 126.04 and 128.56 vs 128.56 — **0.00% error on both**.

Also stamps `var['below_gdna_background']` — the 405 of 5,001 genes (8.1%) whose section-wide
detection rate is at or below background. **Flagged, never dropped**: a gene that is
background-level across the section can still be real in the one cluster that expresses it.

### Stage 2 — `scripts/xenium_celltype/build_phenotyped_base.py`  ✅ IMPLEMENTED, VALIDATED

One QC, one normalisation, one cell population for every downstream stage.

```bash
python3 scripts/xenium_celltype/build_phenotyped_base.py \
    --donors 6533 112 --roi panc --out /data/xenium/processed
```

- Decision #10 unchanged: `MIN_GENES=20`, density quantile 0.98 on `counts/cell_area`
  (never an absolute counts cap — `cell_area` varies 5-10x and a cap removes large healthy
  acinar/ductal cells).
- **`MIN_COUNTS=50` now applies to gDNA-corrected counts.** The number is unchanged; the
  quantity it measures is fixed.
- Preserves the vendor `total_counts` and `genomic_control_counts` in `obs` — these are the
  only record of control counts once the matrix is loaded gene-expression-only, and
  `sc.pp.calculate_qc_metrics` **overwrites** `obs['total_counts']` with a gene-only row sum
  wherever it is called.
- Carries `gdna_lambda`, `gdna_raw`, `gdna_expected_counts` into `obs` so no later stage
  re-opens the bundle.

Measured effect:

| | kept (raw floor) | kept (corrected) | majority-background among kept |
|---|---|---|---|
| 6533 | 96.70% | **91.04%** | 14.84% → **10.75%** |
| 112 | 95.92% | **88.75%** | 18.65% → **13.79%** |

It reduces but does not eliminate background-heavy cells among those kept — a deep cell can be
majority-background and still clear 50 corrected counts. Raising the RAW floor to 300 instead
costs 31% / 34% of cells and still admits them.

### Stage 3 — background-aware detection  ⬜ SPECIFIED

One `detected()` predicate replacing four incompatible definitions. Six decisions in
`analysis/xenium/scripts/hierarchical_phenotype.py` currently rest on **raw count >= 1**:

| constant | value | P(triggered by gDNA alone) |
|---|---|---|
| `IC_COUNT_FLOOR` | 1 | **2.49%** per gene per cell |
| `RAW_DETECT_FLOOR` | 1 | **2.49%** |
| `DOMINANT_CANARY_FLOOR` | 2 | 0.031% — already safe |

**The size bias this creates** — P(a core is spuriously satisfied by gDNA alone, per cell), at
`ANCHOR_MIN_GENES = 2`:

| lineage | core genes | P(>=2 genes by gDNA) |
|---|---:|---:|
| Exocrine | 8 | **0.016** |
| Endocrine | 18 | 0.073 |
| Immune_B | 58 | 0.425 |
| Vascular | 87 | 0.641 |
| Myeloid | 111 | 0.766 |
| Immune_T | 129 | **0.834** |

A 50x differential that scales with core size — large cores are spuriously anchored far more
often. This is a mechanism for the failure `nullgate` measured independently (small cores
Endocrine and Immune_B came out **below** their matched null). Requiring 2 counts per gene
collapses it to 0.0008 for Immune_T.

Fix is primarily a **change of scale, not only a higher floor** — see the lognorm evidence
above. Validate against the Y probes.

### Stage 4 — group-level typing (PRIMARY CALL)  ⬜ reuse `scripts/xenium_celltype/`

`prep -> KODAMA (spatial + nonspatial) -> KMeans -> nullgate`. Already the right architecture:

- **Decides per group, never per cell.** Pooling is the only way to get signal from a
  ~433-transcript cell whose background is ~126 counts.
- **Its null is already detection-matched** (`match_strata`: 100 joint strata on per-gene
  detection rate AND mean). Since gDNA is uniform per probe, a matched random core carries the
  same background as the real core, so it largely cancels in the contrast. Make this explicit:
  record the background-predicted contrast (0 under the null) per group.
- Keep every gate property: `real > max(p95, 0)`, the separability paired null,
  `MIN_GROUP_CELLS = 50`, and refusals that stay refusals.

**`nonspatial` is the only arm that may carry a spatial claim** — a type call that already
knows where the cell is makes "do immune cells sit near islets" partly circular. The arm
difference is the measurement of how much location leaks into identity.

KODAMA is pinned to `smith6jt-cop/KODAMA@d2dbe30e` and `check_shard_equivalence.R` is a release
gate. **Do not substitute a Python embedding** — it discards the equivalence proof.

### Stage 5 — protein arm, reported beside (not merged)  ⬜ SPECIFIED

The same-slide pair gives 70-plex protein on the same cells, far denser than 433 transcripts.
Reuse `core/phenocycler/scripts/celltype/{rules,score}.py` + `phenotype_workflow.csv`
(37 references incl. `Endo-beta/alpha/delta`; 55 of 56 markers map onto this panel).

**Drop the cohort-wide parts** — no UniFORM, no cohort scale. Per-donor floor and robust top,
with the limitation recorded: 6533 and 112 are then not strictly comparable to each other.

The Xenium-only call is primary. **Any protein/RNA concordance claim must come from the
Xenium-only arm**, or it is circular.

### Stage 6 — one decision record, one verifier  ⬜ SPECIFIED

Single schema `{unit, n, label, reason, statistic, null, confidence, runners_up}`. Adopt
`run_donor.composition`'s status vocabulary (`awarded / unavailable / cleared_but_outranked /
no_group_cleared`, **NaN never 0**).

**Fix `verify_run.py` check (8) before any production run.** It tests `lineage ==
"Unavailable"` while `composition` writes `"Unavailable:" + reason`, so it **fails on every
complete run** and `50_assemble.sh` reports "run NOT complete". One-line fix
(`str.startswith`) plus a test.

---

## 3. Islet membership — two backends, one contract

Both donors have islet annotations, but they reached disk differently:

| donor | route | quantity |
|---|---|---|
| 6533 | `Islets_0080870_...geojson` **and** `Parent = Islet_N` in `6533.tsv` | 1,596 Islet + 1,594 IsletExpanded; 95,566 annotated protein cells |
| 112 | `Islets_0081366_...geojson` only (TSV has no islet `Parent`) | 2,378 Islet + 2,376 IsletExpanded |

One `islet_membership(donor)` accessor, two backends, returning
`(cell_id, in_islet, islet_id)` either way. **6533 has both**, so validate polygon containment
against QuPath's own `Parent` on 6533 before applying the geometric rule to 112.

GeoJSON coordinates are **qptiff pixels**; scale is fitted, not read from a tag —
`0.507800 µm/px` with a 24 nm median residual, reproduced independently on both donors.

---

## 4. Traps that have already cost time

1. **Control probes never reach `var`.** `spatialdata_io.xenium` loads Gene Expression only
   (5,001 of 10,029 features). The repo's control-drop regex is a **no-op**, and it would not
   match the genomic controls anyway — they are `Intergenic_Region_*`, and the string
   `Intergenic` appears **zero times** in the monorepo.
2. **`metrics_summary.csv` is read by no code.** Every number in section 0 lives in a file
   nothing opens. Stage 1 now ingests it.
3. **anndata round-trip**: indices must be **unnamed** and every string forced to `object`.
   A named index round-tripped through a second writer reset var_names to '0','1','2'… — the
   identity panel came back 0/5101 and nothing raised. Pinned by
   `tests/test_build_phenotyped_base.py`.
4. **Three mutually exclusive environments**: `xct_py` (no scanpy/scvi/torch), `scvi-env`
   (GPU; `xenium_analysis` is broken for it), `xenium_analysis` + `R_LIBS_USER` for R.
   `conda activate` replaces rather than stacks — the pipeline is a chain of stages, not one
   program.
5. **Don't judge cell typing by proportions.** Repo policy: report raw numbers; no
   "expected fraction" gates. Categorise each cell by what it expresses.
6. **`nullgate_fast.py` (the production path) drifted from the reference once**, silently
   awarding lineages the reference refuses. The frame-for-frame equivalence test must be part
   of the gate, not optional.
7. **`sc.tl.score_genes` already subtracts a binned control-gene reference.** It is the only
   implicit background correction anywhere in the pipeline, it is gene-expression-based rather
   than control-probe-based, and it applies to every `{celltype}_score` column. Do not stack a
   second correction on those scores without checking what it already removed.
8. **Dead config.** `DETECTION_SUBTYPE_LINEAGES[*]["gate_thr"] = log1p(1.0)`
   (`pancreas.py:474, 483, 492, 501`; `thymus.py:724-758`) is **never read** —
   `hierarchical_phenotype.py:870-871` hardcodes `RAW_DETECT_FLOOR` on `layers['counts']`. The
   declared lognorm threshold has no effect, and the comment above it describes behaviour the
   code does not have.
9. **Two independent canary implementations.** `02_phenotyping.ipynb` cell 22 (lognorm,
   `log1p(3)`, exclusive-group counting, `margin >= 2 AND top_depth >= 2`) is a different code
   path from `hierarchical_phenotype.py` (raw counts, floors 1 and 2). The notebook's own
   comment at `:726` says "raw count >= 3" while the code at `:757` is lognorm — and because
   `normalize_total` uses the median library size, a ~40% inflated denominator from gDNA
   raises the effective raw bar. Fix the comment or the code, but know they disagree.
10. **Scripts referenced in `CLAUDE.md` that are missing from disk**:
    `verify_upper_tail_doublets.py`, `lineage_phenotype.py`, `refine_rare_celltypes.py`,
    `compute_panel_detection.py`. `pancreas.py`'s `LABEL_TO_CANARY` is consumed only by the
    first of these, so it is currently dead.

### Highest-leverage locations for the gDNA remediation

`hierarchical_phenotype.py:701-705` (the eligibility mask — one stray transcript makes a
lineage eligible), `:658` (`ic_detect_counts`), `:870-872` (fine gate), `:847-850` (dominant
canary), `01_preprocessing_v2.ipynb:1547` (the vendor `total_counts` overwrite), and `:1657`
(the regex that would need `Intergenic_Region_` / `DeprecatedCodeword_` if the reader ever
stopped filtering controls out).

---

## 5. Validation checklist

| # | check | pass condition | status |
|---|---|---|---|
| 1 | gDNA reproduces the vendor | within 5% of `est_FP..._including_genomic_counts` | ✅ 0.00% both |
| 2 | Y probes under the detection rule | <= ~0.3% (they have no template) | ⬜ Stage 3 |
| 3 | core size bias | P(Immune_T satisfied by gDNA) 0.834 -> <0.01, Exocrine ~unchanged | ⬜ Stage 3 |
| 4 | null calibration | permuted labels award at the pre-registered rate (1 lineage 6.2%, 8 lineages 30.2% @ 50 draws) | ⬜ Stage 4 |
| 5 | islet containment | geometric rule reproduces QuPath `Parent` on 6533 | ⬜ |
| 6 | arm difference | `spatial - nonspatial` per k and lineage, never averaged | ⬜ |
| 7 | `verify_run.py` | passes on a complete run (currently cannot) | ⬜ |

---

## 6. Reference numbers

```
panel            5,001 genes | 21 genomic_control | 40 negative_control | 3,291 deprecated
                 genomic controls are Intergenic_Region_*

0080870 (6533)   717,461 cells | median 433 transcripts | lambda 0.02520 | 126.04 FP/cell
0081366 (112)    671,211 cells | median 402 transcripts | lambda 0.02571 | 128.56 FP/cell

P(>=1 count from gDNA alone)   2.49% / 2.54% of gene-cell pairs
P(>=2 counts)                  0.031%
non-gDNA floor (Y probes)      ~0.3%
genes at/below background      405 of 5,001 (8.1%)
cells majority-background      15.53% / 19.77%
cells where bg > all counts     2.82% /  3.51%
px -> um scale (PhenoCycler)   0.507800, residual 24 nm (fitted, both donors)
```
