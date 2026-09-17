# Handoff — rebuild the PhenoCycler typing comparison on a sound substrate

**For a new session.** Written 2026-09-08 against `main` @ `b5e2ecf`. You are picking up work whose previous
conclusion has been **withdrawn**. Read this document, then the four files in "Start here", before running
anything. Nothing is currently running; no process needs babysitting.

---

## 1. Start here

| read | why |
|---|---|
| `docs/BROADTYPE_PROTOCOL_2026-09.md` §4.9 | the withdrawal and the twelve verified defects. This is the ground truth for what went wrong |
| `core/phenocycler/scripts/broadtype/gate.py` | the retired decision logic — you need to recognise it to avoid rebuilding it |
| `core/phenocycler/scripts/celltype/README.md` §alignment | the repo's own rule that six identity markers may not take a cohort-wide cutoff |
| `core/phenocycler/DATA_README.md:243` | the repo's own rule that `analysis_eligible` gates the identity call |

Repo rules that still bind and are easy to trip over: `.venv` is never recreated (a Globus daemon runs under
it); `git reset --hard` and `git clean` are unsafe (~80 GB of untracked ignored data); nothing writes
`data/runs/`; editing `core/phenocycler/phenocycler/**` re-keys a future run; `/home` is 86% full so pixel
scratch goes to `/data`; QuPath classifications are audit context only.

## 2. State of the work

**Withdrawn.** The producer choice. `data/development/broadtype/release/` still exists on disk and must not be
read as a verified typing. The registry latch `NO_VERIFIED_TYPING` never moved, so nothing downstream was
released on it.

**Stands and is reusable.** The QC figures (P01–P05, X01–X03) are measurements of the run. The integration
pipeline ran all fourteen stages for both donors and its outputs are on disk. The Xenium release exists. The
report lane, its manifest and `--verify` work (62 figures, 261 files, no drift). The mechanical bookkeeping —
`verdict.py` recomputing survival from records, `gates.require_evidence` failing closed, the tiering rule —
is sound and worth keeping.

**Superseded by this handoff.** `docs/BROADTYPE_PROTOCOL_2026-09.md` Section 2 (the sampled frame) and
Section 3 (single-marker removal, `tail_q`, `heldout_frac`, `subsample_frac`). The protocol document itself
stays as the record; do not edit above its results line.

**Two representative donors, unchanged:** 6436 (T1D) and 6442 (ND).

## 3. The five rulings this handoff implements

Given by the maintainer on 2026-09-08. They are not open for re-litigation; if one turns out to be
unimplementable, stop and say so rather than working around it.

1. **All cells, every donor, every method.** Pancreas tissue is not homogeneous, so a sampled frame is not
   evidence about a donor.
2. **No cell subsampling anywhere in validation.** Heterogeneity confounds a subsample used to check quality.
3. **Single-marker removal is out.** Panel design defines a broad type by a combination of markers.
4. **The compute estimates were wrong, the hardware is under-used, or both.** Treat the measured runtimes as
   evidence of bad execution, not of intrinsic cost. No method is dropped for being slow.
5. **Exclusions have reasons.** Find and report them; do not reverse them silently.

### How each ruling changed the plan

**Ruling 1 and 2 together retire far more than the frame.** The old design sampled twice over: sixteen
windows out of a slide, and then, inside those windows, statistics computed on subsamples — a 42-cell marker
tail, a 200,000-cell cap on the Moran graph, and 90% subsample replicates for stability. All of it goes.
Phase 3 is rebuilt from nothing rather than rescaled.

They also force a change nobody would notice otherwise. At whole-donor scope the existing spatial null
permutes labels *globally*, which destroys tissue architecture and collapses the null p95 from ~0.005–0.029
to ~0.0014 for every label — the coherence test becomes vacuous exactly when you scale it up. Phase 3
therefore specifies a null that permutes within 500 µm blocks.

**Ruling 3 removes the blocking test and, with it, most of the compute.** Twenty-five of thirty-one runs per
donor per method existed only to remove one marker at a time. Phase 3 replaces them with a criterion built on
the combination the panel already declares, and Phase 2 inherits a 5× smaller workload.

**Ruling 4 sent me back to the measurements rather than to a scoping decision.** I had proposed dropping
methods or accepting four days of GPU. The audit found four specific causes instead, all of them fixable, and
Phase 2 is now about execution. No method is dropped.

**Ruling 5 removed a proposal I should not have made.** I had offered to recompute the spillover correction
keeping the excluded polygons. The exclusions are principled — see §4 — so the plan now applies the repo's own
rule consistently and reports the consequences instead of reversing them.

**A separate correction.** I had called FOXP3 and TP63 anomalies. They are not: all eight
nucleus-compartment markers are handled identically, and those two are simply the only ones that are identity
markers. The real defect is that the identity vector mixes compartments and compensation states, which is what
Phase 1 fixes.

## 4. The exclusion reasons, as asked

| donor | flag | excluded | dominant reason | remainder |
|---|---|---|---|---|
| 6436 | `analysis_eligible` | 353,357 (18.3%) | `missing_nucleus_geometry` 249,390 | `cell_area_below_min` 109,387; `low_solidity` 14,896; `cell_area_above_max` 665; 22 raster |
| 6436 | `spillover_context_eligible` | 249,407 (12.9%) | `missing_nucleus_geometry` 249,390 | 21 raster |
| 6442 | `analysis_eligible` | 295,614 (17.2%) | `missing_nucleus_geometry` 225,325 | `cell_area_below_min` 70,966; `low_solidity` 15,804; `cell_area_above_max` 555; 24 raster |
| 6442 | `spillover_context_eligible` | 225,345 (13.1%) | `missing_nucleus_geometry` 225,325 | 24 raster |

A detection with no nucleus polygon is the whole story: 71% and 76% of the analysis exclusions, 99.99% of the
spillover exclusions. Such an object has no nucleus in which to measure the eight nucleus-compartment
markers, and `redsea.apply_spillover_context_eligibility` (`phenocycler/redsea.py:661-678`) erases its pixels
from the label image before any channel is read, so every corrected column is NaN. Both are coherent.

Two consequences are **reported, not fixed**, because fixing them means re-streaming 82 GB qptiffs with no
saved intermediates. Erased polygons are absent from the contact graph, so surviving neighbours keep the
boundary reinforcement and lose the matching subtraction — measured on 6436, the corrected/raw ratio runs
0.92 with no erased neighbour to 1.05 with three or more. And the erasure is spatially uneven, 9.8% to 25.8%
across the old frame regions.

---

## Phase 0 — freeze a new protocol before anything runs (0.5 day)

`docs/TYPING_PROTOCOL_2026-09-08.md`, sections 0–3 committed **before the first record**, same append-only
idiom as its predecessor, which it cites and supersedes in part. New lane
`core/phenocycler/scripts/typing/`; say *method*, not *arm*. Criteria live in a new `criteria.py` carrying
`assert_frequency_free` and a new `assert_no_cell_subsampling` that fails on any key naming a fraction or cap
of cells.

## Phase 1 — one substrate every method reads (1.5 days)

`scripts/typing/substrate.py` → `data/development/typing/substrate/donor_id=<d>/data.parquet`, one row per
QuPath object in expression order, written through `phenocycler.expression.write_donor_partition` (atomic,
refuses overwrite) and fingerprinted with `artifacts.StageManifest`. The current uniform tables are ragged
(57 vs 58 markers between donors) and unfingerprinted — nothing on disk can prove which bytes a method read.

- All 28 identity markers on **one measurement basis**: FOXP3 and TP63 move from `cells.Nucleus__<m>` to
  `redsea.Cell__<m>`, which exists for all 20 donors, so every identity column is a whole-cell mean. This
  fixes the compartment mix and not the compensation one — `Cell__FOXP3` is passthrough, 1.8% from the QuPath
  whole-cell mean — and that residual is recorded on the marker's row rather than hidden.
- `eligible_for_identity` = `analysis_eligible`, applied **once, here**. Ineligible cells stay in the table
  and type as `Unavailable(<first reason>)`, which is what `DATA_README.md:243` prescribes and what
  `release.build_donor_calls:362-373` already does at the wrong end of the pipeline.
- NaN stays NaN. Replace `standardise_within_donor`'s `np.where(isfinite, Z, 0.0)`, which silently places
  ~13% of cells at the marker mean.
- A declared transform applied here rather than left to each method: `arcsinh(x / cofactor)` from each
  marker's own measured background, then within-donor standardisation. On the current linear scale the
  z-score IQR spans 41.6× across markers and runs *inversely* with specificity — SMA contributes 0.17% of a
  typical Euclidean distance while CD79a contributes most of it.
- The six identity markers whose alignment failed its own 0.10 bar (Vimentin 0.852, Pan_Cytokeratin 0.405,
  E_cadherin 0.323, CD31 0.155, INS 0.145, Ker8_18 0.109) carry an `over_bar` flag and a within-donor rank
  column beside the aligned value. `scripts/celltype/README.md:249` already rules that a single cohort-wide
  cutoff is not valid for them; any method using one must read the rank.
- **Pixel methods read the same cell set**: also emit `masks/donor_id=<d>.npy`, the whole-slide `int32` label
  raster with ineligible labels zeroed, so an image method physically cannot see a cell the table methods
  refuse.

Reuse `scripts/kintsugi/features.py:184-211 load_block` as the reference filter and
`artifacts.snapshot_parquet_dataset` as the validator.

## Phase 2 — whole-donor execution, with the machine actually used (2 days)

All 1,932,797 / 1,715,283 cells, every method. `--scope donor` already works for the seven table methods and
two have already run there. The four pixel methods need whole-slide tiling, and their measured cost is an
artefact of how they were driven:

| what was done | cost | what to do instead |
|---|---|---|
| 2048-px windows around 1024-px cores | every cell inferred ~4×; 179,770 inferred to score 48,064 | non-overlapping 2048-px tiles over the 771 / 668 populated tiles; ~40% of the grid is empty and skippable |
| KRONOS2 at batch 16, GPU at 3–10% of a 19.3 GB card | 46–106 ms/cell | raise the batch until the card is full; batch 16 is for bit-exactness, not required here; record what was used |
| embedding cache keyed on the channel set | each of 25 marker removals re-embedded the whole donor | ruling 3 retires them; 25 of 31 runs per donor disappear |
| mask read from Deflate TIFF, RowsPerStrip=2 | 244 MB decoded to keep 16 MB | hold the whole-slide raster in RAM: 5.66 / 4.64 GB against 490 GB available |

Whole-slide rasterisation is measured at 163 s / 155 s (`redsea.rasterize_mask`); re-run
`window_mask.verify_against_archive` over the **full** slide before trusting the archived masks, which are
currently proven over 67 Mpx of 1,414. Reusable: `common/tiff_window.read_window`,
`restain/tiles.read_window_on_grid`, `restain/tiles._sparse_chunk_mask`, `restain/explorer._tile_stream`. One
`TiffFile` handle per worker process. Scratch to `/data`.

Three methods were also driven in configurations their authors would reject, and that is corrected here:
deepcell-types **not** on the channel set this repo already measured as its worst
(`docs/DEEPCELL_TYPES_REVIEW_2026-08-24.md`); the SOM **with** KINTSUGI's meta-clustering step
(`som_adapter.meta_cluster` exists and was never called); Nimbus with an input normalisation that does not
call FOXP3 positive in 68% of cells; TACIT on a normalised scale, since its internal `slope < 0.05` demotion
constant assumes one.

## Phase 3 — validation with no subsampling and no single-marker removal (2 days)

`scripts/typing/gate.py`. Three blocking criteria, each on **every** cell of the donor.

1. **Combination agreement.** The declared broad type already *is* a combination: `rules.requirements`
   (`scripts/celltype/rules.py:290`) returns `(all_of, [any_of groups], none_of)` per label and
   `score.meets_requirements` (`score.py:180`) evaluates it per cell against each marker's own measured
   background, failing closed on NaN. Both exist and both already run on full donors. The criterion is the
   agreement between a method's label and the declared combination over all cells, reported as a full
   confusion matrix, against a null that permutes the method's labels across all cells.
2. **Spatial coherence, all cells, blocked null.** Moran's I of the one-hot label on the full k = 12 graph
   over every eligible cell, no cap; null permutes labels within 500 µm blocks, not globally.
3. **Repeatability.** Stochastic methods re-run at different seeds on the full cell set. No 90% draws.

Retired and stated as retired: `tail_q`, `tail_min_k`, `min_tail_in_pool`, `heldout_frac`, `subsample_frac`,
`spatial_max_cells`. `Ambiguous` and `Unavailable` are counted identically by every statistic, so a method
cannot change its verdict by relabelling abstentions. Tests are corrected for the number of labels offered,
and the count performed is printed beside the table.

## Phase 4 — bidirectional cross-modal evidence (1.5 days)

The independent check the old design lacked entirely — `scripts/broadtype/gate.py` contains no reference to
Xenium. Both directions, every donor, all cells.

- **PhenoCycler → Xenium.** Composition over the seven comparable lineages (`integration/donor.py
  compare_donor`), and per-anchor protein↔gene concordance. The latter needs `donor.py:164` fixed: it reads
  `feat__score_<lineage>`, constant within a lineage, so the statistic is 52 copies of one number. The 26
  per-gene `feat__<GENE>` columns it should use are already on disk in `cells_xen`.
- **Xenium → PhenoCycler.** The Xenium tier-1 rule is itself a combination test — a clean anchor detects ≥2
  genes of exactly one identity core — which makes it the natural independent reference. Compare in the
  shared 7-label vocabulary that `vocab.XENIUM_TO_COMMON` and both release label spaces already agree on.
- **Registration-free niche agreement, both directions**, via `grid.niche_profiles` / `niche_agreement` with
  the 1000-draw lineage-permutation null already implemented in `report/multi.py:m06_niches`. This is the
  strongest cross-modal number in the repo today (per-niche ρ 0.829 / 0.737 against p95 0.46 / 0.44) and it
  works at any registration verdict, which matters because 18 of 20 donors have no registration.
- **Same-section IF**, the only place a protein and a transcript are measured on the same cell:
  `feat__if_INS` / `feat__if_GCG` on 100% of Xenium cells behind a registration that passes its own bar
  (0.55 / 0.63 µm median). It speaks to **Endocrine only** — no other panel protein has been read on Xenium
  cells, Opal 570 failed on every nPOD-batch section — and 6436's INS channel is unusable (2.2× islet
  enrichment against GCG's 66×), so GCG carries 6436 and INS carries 6442. Say so on the figure; never
  average the donors.
- Registered statistics run only where `qc.ALLOWED_BY_VERDICT` permits: 6442 WARN, 6436 FAIL. Report both.

Cross-modal agreement blocks **alongside** the Phase 3 criteria. A label verified within the panel but
contradicted across platforms is reported as contradicted, not as verified.

## Phase 5 — records and figures (1 day)

Append results below the new protocol's results line. Re-render through
`scripts/report/run_report.py --section all --verify`; keep the withdrawn figures but mark them superseded.
Update `docs/WORKFLOW_2026-09-02_MAP.md`, `docs/EXTERNAL_PACKAGES.md` and the memory note. Write a release
**only** if a method survives; otherwise record "no producer chosen" and leave the latch down.

## Verification

- `pytest -q scripts/typing/tests`, plus existing `scripts/report/tests`, `tests/integration` (308),
  `python -m pytest -q`.
- New tests pinning the rulings: `test_no_cell_subsampling_in_criteria`;
  `test_every_method_sees_one_cell_set` (record-level: identical eligible cell count across methods per
  donor); `test_identity_columns_are_one_measurement_basis`; `test_spatial_null_is_blocked_not_global`.
- Substrate: `snapshot_parquet_dataset` passes; rows 1,932,797 / 1,715,283; NaN count equals the
  spillover-context exclusion exactly; `eligible_for_identity` sums to 1,579,440 / 1,419,669.
- Masks: `verify_against_archive` reports 0 differing pixels over the full slide.
- Execution: every record shows `scope = donor` and `n_cells` equal to the donor's object count.
- Cross-modal: `marker_gene_concordance` returns a per-anchor statistic that varies within a lineage.

## Not in this handoff

Reversing any geometry-QC exclusion. Re-running REDSEA to retain erased polygons (~4.9 M cells cohort-wide,
no saved intermediates, a full re-stream of every qptiff). Re-keying the run. Renaming the `arm` identifiers
already on disk. Expanding beyond the two donors. KRONOS v1, whose Hugging Face access was never granted.

## Traps that cost time in the previous session

- A `pgrep` watcher whose own command line matches its pattern never exits. Three shells looped for up to two
  days. Match on a pid, or use the harness's own waiting.
- `data/development/broadtype/` is 24 GB and `/home` is 86% full. Pixel scratch and embedding caches must go
  to `/data`.
- The figure manifest used to keep the first hash it saw for an input, so `--verify` reported drift on stale
  records. It now re-fingerprints on every read; keep that.
