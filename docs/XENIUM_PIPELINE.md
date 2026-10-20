# Xenium analysis pipeline — gDNA-corrected, per-donor

**Status:** Stages 1, 2 and 5 implemented and validated. Stages 3, 4 and 6 specified, not yet built.
**Order:** `1 gdna -> 2 QC -> 5 protein arm -> 3 context-conditional detection -> 4 RNA group
typing -> 6 records`. The protein arm runs **before** the RNA arm, per the maintainer: it is
the dense measurement (70-plex on the same cells against ~433 transcripts of which ~126 are
background), so it is the natural source of the context prior the RNA denoising needs, and
running it first means that prior does not depend on the RNA.
**Scope:** Xenium Prime 5k, per donor. Nothing in this pipeline is cohort-wide.
**Applies to:** the Pheno-Xeno pair — 0080870 (donor 6533) and 0081366 (donor 112).
Both donors are **female** (confirmed by the maintainer).

> **Standing note, 2026-09-18 — what happened to these two donors while stages 3-6 stayed
> unbuilt.** Both are typed in the 21-donor pancreas release
> `xenium-panc-cluster-name-20260917` (`docs/PANC_XENIUM_TYPING_2026-09-17.md`), produced by a
> different method — per-section Leiden, naming-only core scoring, `lognorm >= log1p(1)` detection,
> no raw-count floor — and **with no gDNA correction applied and no flag on the rows**. The
> correction in §1 was never run over that release. Any call on 6533 or 112 taken from it, and any
> cohort-pooled number that includes them, carries the §0 background silently. The cohort-wide
> triage that partitions the rest of the cohort against this pair is
> `docs/XENIUM_GDNA_TRIAGE_2026-09-17.md`; the pLN half, measured on λ rather than the adjusted
> rate, is `docs/PLN_IMMUNE_TYPING_HAZARD_2026-09-17.md`.

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
cd core/multiplex5000
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
  exocrine cells). *(Corrected 2026-09-18: this read "acinar/ductal cells". The 5K panel carries
  no PRSS1, CPA1, CELA3A, CTRB1, CTRB2, AMY2A or REG1A — measured per section — so **Exocrine is
  one label and no acinar or ductal claim is made or supported here either**. The QC rationale is
  unchanged; only the label is.)*
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

KODAMA is pinned to `smith6jt-cop/KODAMA@d2dbe30e` and `check_shard_equivalence.R` was the release
gate (quarantined 2026-09-23 under data/ARCHIVE_2026-09-15/quarantine/ with the rest of `xenium_celltype/hpg/`). **Do not substitute a Python embedding** — it discards the equivalence proof.

### Stage 5 — protein arm  ✅ IMPLEMENTED, VALIDATED

`core/multiplex5000/scripts/celltype/run_sameslide_protein.py`. Reads the same-slide lane's own
`expression/donor_id=<d>/data.parquet` — REDSEA spillover-corrected whole-cell means for the
corrected markers, raw QuPath *nucleus* mean for the eight in `panel.NUCLEUS_MEASURED`,
verified column by column against `redsea/` and `cells/`.

```bash
python core/multiplex5000/scripts/celltype/run_sameslide_protein.py            # both donors
SAMESLIDE_PHENO_DATA=<dir> python .../run_sameslide_protein.py --donors 6533
```

37 references over 56 rule markers, `signed` encoding, `require_min = 0.10`. Output:
`data/development/celltype/sameslide_protein/donor_id=<d>/{data.parquet,scale.csv}`.
Runtime ~8 s per donor.

**Cohort-wide parts dropped.** No UniFORM (the Step-2 aligned partitions exist only for the
22 cohort donors, and these two are not among them) and no cohort scale. Whole-cell INS means
differ 11.6 (6533) against 113.0 (0112), so one pooled floor would put an entire donor under
it. Per-donor floor (mode of `log1p`, the same estimator `build_scale` uses) and robust top
(99.9th percentile). **The cost is real and is recorded in the module: within a donor the
scaled values are comparable, between donors they are not.**

**`Ker8_18` is on the rule table and not on this panel**, emitted as an all-NaN column.
`audit_unmeasured` separates the three consequences rather than conflating them, and prints
all three per donor:

| consequence | references | verdict |
|---|---:|---|
| OR-group loses one of four alternatives (`E-cadherin \| Ker8/18 \| Pan-CK \| EpCAM`) | 11 epithelial | gate still satisfiable |
| negative constraint lost | 16 immune | less specific, reported |
| mandatory constraint or whole OR-group lost | 0 | none dropped |

**Inside `Islet_N` the decision is an ordered cascade, not a competition.** The polygons were
drawn by a QuPath pixel classifier on INS, GCG, SST and CD99, so the cascade reads those four
in the order the annotation was built:

| step | test | result |
|---|---|---|
| 1 | hormone argmax over INS/GCG/SST, clearing the floor | `Endo-beta` / `Endo-alpha` / `Endo-delta` |
| 2 | CD99 positive, for cells with no hormone | `Endo-hormone-negative` |
| 3 | centroid inside a drawn structure **and** its marker | `Endo-vascular` / `Endo-lymphatic` / `Neural` |
| 4 | marker profile, non-endocrine references | everything else |

**Shape** is `geometry_qc`'s `min_cell_solidity = 0.70` (`config.ini:60`), applied upstream as
part of `qc_analysis_eligible` — 4,144 (6533) / 2,169 (0112) cells carry the `low_solidity`
reason. **Structure** is read from the QuPath object layers (`Blood_vessel`, `Lymphatic`,
`Nerve` — 111,806 / 8,535 / 11,293 polygons on 6533), point-in-polygon exactly as islet
membership is. Step 3 requires **both**: a cell inside a vessel showing no CD31 is not a
vessel cell and falls through to step 4. On 6533's cores that is 534 of the 2,901 cells inside
a vessel polygon; on 0112, 1,089 of 4,217.

The layers and the marker profile are independent and agree far above chance, measured
**outside** the islet where the cascade never consults a layer:

| layer | called its reference | outside any layer | fold |
|---|---:|---:|---:|
| `Blood_vessel` | 37.5% / 29.9% | 0.9% / 0.6% | 41x / 47x |
| `Lymphatic` | 36.4% / 10.6% | 0.04% / 0.01% | 921x / 1177x |
| `Nerve` | 50.9% / 32.1% | 0.5% / 0.2% | 107x / 189x |

Step 4 scores the non-endocrine references from `run_islet.intruder_references`, which excludes
the whole epithelial branch, ancestors and siblings alike — an endocrine cell satisfies those
rules by construction. A hand-rolled filter excluding only endocrine-prefixed names was tried
first: the generic `Epithelial` parent then took 21,450 of 65,363 islet cells.

**The lane's step order holds.** Step 1 (marker usability, per donor) runs first and gates the
markers scored here; `run_step1.py` was made runnable on a donor outside the acquisition table
(`CELLTYPE_EXPRESSION_DIR`, plus on-disk marker resolution) rather than substituted for, and it
rejects **0 of 58 markers** on both — the same answer it gave on all 1,120 cohort (donor,
marker) cells. Step 2 is the only step that cannot run, and not for convenience: it is UniFORM
alignment **across the cohort**.

**Measured — core cells, median scaled hormone by call:**

| call | INS | GCG | SST | n (6533 / 0112) |
|---|---:|---:|---:|---|
| `Endo-beta` | **0.637** | 0.082 | 0.055 | 2,319 / 6,350 |
| `Endo-alpha` | 0.000 | **0.405** | 0.054 | 6,505 / 3,024 |
| `Endo-delta` | 0.000 | 0.095 | **0.524** | 2,250 / 2,001 |
| `Endo-hormone-negative` | 0.000 | 0.031 | 0.019 | 1,247 / 3,655 |

Each endocrine call carries its own hormone as its median and not the others. Fraction of
endocrine-by-prior cells clearing the 0.10 gate: **89.9% (6533), 75.7% (0112)** — against
**37.9%** on uncompensated intensities, which is what REDSEA was being waited on for.

**Fixed, and what it cost.** An earlier arrangement removed the endocrine references from
contention inside the core and made a cell endocrine only by meeting NO other rule. It had no
endocrine escape hatch and the measurement said so: 98.9% of the cells it called
`Imm-NK-cell` inside 6533's cores cleared the hormone floor, 1,969 of them SST-top, at a
median Granzyme B and CD3e of 0.00. A delta cell satisfying `CD56` + (`CD57`|`Granzyme_B`) +
any one of nine pan-immune markers had nowhere else to go. Under the cascade:

| 6533 core (n=26,388) | before | after |
|---|---:|---:|
| `Imm-NK-cell` | 3,002 | **6** |
| islet-endocrine | 12,321 | **24,064** (91.2%) |
| `Endo-alpha` / `delta` / `beta` | 6,505 / 2,250 / 2,319 | 11,794 / 7,177 / 3,532 |

0112: `Imm-NK-cell` 1,814 → 43; islet-endocrine 15,030 → 22,774 (73.1%). Every non-endocrine
call left in a core sits at a median scaled CD99 of 0.000–0.015 and a median top hormone of
0.00–0.04, so step 4 operates on cells that are genuinely neither.

The Xenium-only call is primary. **Any protein/RNA concordance claim must come from the
Xenium-only arm**, or it is circular.

### Stage 6 — one decision record, one verifier  ⬜ SPECIFIED

Single schema `{unit, n, label, reason, statistic, null, confidence, runners_up}`. Adopt
`run_donor.composition`'s status vocabulary (`awarded / unavailable / cleared_but_outranked /
no_group_cleared`, **NaN never 0**).

**Fixed 2026-09-17 (990c5b2), kept as a record.** `verify_run.py` check (8) tested
`lineage == "Unavailable"` while `composition` writes `"Unavailable:" + reason`, so it failed
on every complete run and `50_assemble.sh` reported "run NOT complete". It now uses
`str.startswith` (`verify_run.py:163`), with `tests/test_verify_run.py`. The `hpg/` run
apparatus it verified, `50_assemble.sh` and the rest, was quarantined 2026-09-23 under
data/ARCHIVE_2026-09-15/quarantine/.

---

## 3. Islet membership — two backends, one contract

Both donors have islet annotations, but they reached disk differently:

| donor | route | quantity |
|---|---|---|
| 6533 | `Islets_0080870_...geojson` **and** `Parent = Islet_N` in `6533.tsv` | 1,596 Islet + 1,594 IsletExpanded; 95,566 annotated protein cells |
| 112 | `Islets_0081366_...geojson` only (TSV has no islet `Parent`) | 2,378 Islet + 2,376 IsletExpanded |

One `islet_membership(donor)` accessor, two backends, returning
`(cell_id, in_islet, islet_id)` either way — implemented as
`run_sameslide_protein.islet_regions`, which uses the `Parent`-derived `cell_region` when the
partition carries one and recomputes containment when it does not.

**The geometric route is validated against QuPath's own `Parent` on 6533**, the donor that has
both, which is the gate this section asked for before applying it to 112:

| | QuPath `Parent` | geometric | Jaccard |
|---|---:|---:|---:|
| core | 30,197 | 31,001 | 0.926 |
| peri | 65,308 | 64,486 | 0.922 |
| islet-associated (core \| peri) | 95,505 | 95,487 | **0.961** |

**99.44% exact agreement over all 1,027,546 cells.** The disagreement is a boundary effect —
1,358 QuPath-peri read geometric-core and 569 the other way, against 206 + 221 that change
islet-associated status at all. 112's own core:ring ratio from the geometric route (35,656 :
81,517 = 2.29) lands beside 6533's QuPath ratio (30,197 : 65,308 = 2.16).

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
| 3b | geometric islet route vs QuPath `Parent` (6533) | Jaccard on islet-associated | ✅ 0.961, 99.44% exact |
| 3c | endocrine call carries its own hormone | each of INS/GCG/SST is the median of exactly one call | ✅ Stage 5 |
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
