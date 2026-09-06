---
title: Broad-typing bake-off protocol (PhenoCycler)
status: pre-registered 2026-09-05; Sections 0-3 are frozen before the first screen record
supersedes: none (the DeepCellTypes review of 2026-08-24 stays closed and is cited, not reopened)
creates_cell_calls: development only (data/development/broadtype/); no production output
production_outputs_written: none; registry latch NO_VERIFIED_TYPING unchanged
owner: smith6jt
lane: core/phenocycler/scripts/broadtype/
---

# 0. What is decided here, before any result

The maintainer's request (2026-09-04): systematically test and compare the cell-typing packages in the
`smith6jt-cop` GitHub, including KRONOS, on the BROAD cell types the panel's best markers can separate;
then finalize. "Finalize" means a written producer choice (plan Track A.3) and a development-stamped
release contract (A.4a). The registry latch stays down.

Everything in Sections 1-3 is written before the first arm runs. Results are appended below the line
in Section 4 and never edited into the bar. Any configuration change to an arm after its result is
seen is a NEW arm name. If no arm survives, no producer is chosen and the comparison table ships.

Rules in force: root `CLAUDE.md`; no frequency / uniformity / disease-ordering tuning; QuPath classes
and legacy labels are audit context only; `object_id` is the join key; every conclusion is scoped to
"this arm at this configuration on donors 6436 and 6442".

# 1. Label space and tiering

`scripts/broadtype/declarations.py` derives the label space by code from two pre-existing
declarations -- the marker registry (`kind`, `broad_compartment`) and the 37-reference rule table
(only for the Immune branch split and the endocrine hormones) -- and writes `declarations.json`
plus `rules_broad8.csv`; `tests/test_declarations.py` regenerates both and fails on a hand edit.

Tiering rule (the PhenoCycler mirror of Xenium decision #17): identity inputs are registry `type`
markers only; `state` markers (23) never enter identity and are appended after typing as continuous
`state_<marker>_donor_rank` columns; nothing is thresholded into a positive call. Conflict rule: a
marker whose rule-table use makes a definitional (`pos`/`allpos`) identity claim in a different
top-level branch than its registry compartment is excluded; OR-membership elsewhere is a note.

Derived on 2026-09-05 (registry sha and rule-table sha recorded in `declarations.json`):

| lineage_native | common | declared identity markers | held-out gated |
|---|---|---|---|
| Endocrine | Endocrine | GCG, INS, SST | no -- mutually exclusive subtypes; gated by spatial coherence + stability |
| Epithelial_nonendocrine | Exocrine | E_cadherin, EpCAM, Ker8_18, Keratin_5, Pan_Cytokeratin, TP63 | yes |
| Vascular | Vascular | CD31, Caveolin, Podoplanin | yes |
| Stromal | Stromal | SMA, Vimentin | yes |
| T_NK | T_NK | CD3e, CD4, CD8, FOXP3 | yes |
| B_Plasma | B_Plasma | CD20, CD79a | yes |
| Myeloid | Myeloid | CD11c, CD163, CD206, CD209, CD68, Iba1, M2Gal3, MPO | yes |

Excluded by the conflict rule, with the row that caused it: B3TUBB (registry Neural; `allpos` in
`Endo-hormone-negative`), CD34 (registry Endothelial; `allpos` in `Mes-CD34-fibroblast`), CD66
(registry Epithelial; `allpos` in `Imm-neutrophil`), CD56 and b_Catenin1 (no registry compartment),
CD11b and IAPP (registry `excluded`). **Neural is therefore not declared**: its only registered marker
fails the conflict rule. That is a consequence of the pre-registered rule, not of any result, and the
protocol keeps it. The plan of 2026-09-05 listed Neural as "reported, not gated"; this document is the
correction.

No output of this lane says acinar or ductal: the panel declares no acinar- or ductal-specific marker.

# 2. Frame and cell set

`scripts/broadtype/frame.py`: 16 hash-drawn 1024-px cores per donor (8 blind + 8 islet-anchored) by the
DCT lane's `regions.select_regions` rule with `max_window = 2048`; the islet draw uses QuPath parentage
as geometry only. One mask-derived cell set per donor for every arm: labels inside each core from
`window_mask.rasterize_windows` (label = feature index + 1, so `object_ids[label - 1]` is the UUID),
bit-identity against the archived whole-slide masks (`verify_against_archive`) as a stop-the-lane
gate. `tails_in_frame.csv` records, before any arm runs, how many Step-1 tail cells each marker has in
the frame; a marker below `min_tail_in_pool` is judged whole-donor instead. Table arms run whole-donor
as well; image arms run on the frame and whole-slide only as survivors.

# 3. Arms, criteria, verdict

Arms (each `scripts/broadtype/arms/<name>.py`, one common record schema):
A `cosine-broad8`, A' `cosine-37-collapsed` (continuity row; inherits the state-gated table, flagged),
B `som-kintsugi` (type-marker block), C `tacit`, D `celesta` (default thresholds; its equal-proportions
mixture prior flagged), E `nimbus-rules`, F `kronos-clusters` (conditional on Hugging Face access; the
request is pending review as of 2026-09-05), G `dct-ptft-broad8-typeonly`, H `tasp` (tail-anchored kNN
propagation), I `consensus` (descriptive only). Every arm's identity input is type markers only;
`test_no_state_marker_in_identity_inputs` asserts it from the records.

`scripts/broadtype/criteria.py` is the bar, verbatim: tail_q 0.001, tail_min_k 20, min_calls 50,
min_tail_in_pool 20, null_draws 200, null_quantile 0.95, morans_k 12, n_replicates 5, primary_seed 0,
subsample_frac 0.9, n_blind 8, n_islet 8, core_px 1024, window_px 2048, heldout_frac 2/3,
survive_frac 2/3. `assert_frequency_free` runs at import.

Per arm, per donor (pooled over the frame with within-region nulls, or whole-donor):
1. Held-out marker enrichment (blocking; labels with >= 2 declared markers): withhold marker m from
   the arm's INPUT, re-run; h = P(label | m-tail cell) / P(label | judged cell) against 200
   size- and region-matched random cell sets; clear iff h > null p95 with >= 20 tail cells and >= 50
   calls. Label verified iff >= ceil(2/3 x n_judged) markers clear and n_judged >= 2.
2. Spatial coherence (blocking; every label): Moran's I of the one-hot label (k = 12) against 200
   within-region label permutations; clear iff I > p95 with >= 50 calls.
3. Stability (blocking; stochastic arms): the spatial verdict is invariant across 5 seeds or 5 x 90 %
   subsample draws; ARI is descriptive.
4. Cross-lineage tail contamination, QuPath islet-core enrichment of Endocrine (AUDIT ONLY), cross-arm
   agreement, and -- once integration has run -- cross-modal niche/islet concordance: descriptive.
Unjudgeable is never a pass. Null-gated cluster annotation (arms B, F): a cluster is awarded a label
iff its declared-marker contrast beats the p95 of 200 size-matched random cell sets AND ranks first by
a 25 % margin; the false-award rate on synthetic noise is measured by the tests and printed beside
every table.

Verdict: an arm survives iff >= ceil(2/3 x 7) = 5 labels are verified on BOTH donors. Producer
choice (A.3), mechanical: among survivors, (i) verified-label count, (ii) mean held-out margin over
the null p95, (iii) fewest external pins / lowest whole-donor runtime. Agreement with another arm or
modality never chooses an arm, a k, or a threshold.

---

# 4. Results (appended below this line, never edited into the bar)

## 4.0 Arms declared after the first records (2026-09-06)

Per Section 1.3, a configuration change after results is a NEW arm name; the original rows stand.

- **`som-kintsugi-aligned`** — the plan's literal arm-B specification (KINTSUGI `features.py` arm-A block:
  aligned intensities, no rank; StandardScaler inside `som_cluster`; type markers only). Declared after
  `som-kintsugi` (within-donor-rank features) had been screened on both donors and awarded no lineage
  cluster on either (its clusters separate by overall brightness). Everything else identical to arm B.
- **`kronos2-coral`** — KRONOS2 (marker-aware DINOv2 ViT-B/16, Hugging Face `MahmoodLab/KRONOS2`
  snapshot `e2a54732`) per-cell CLS embeddings obtained THROUGH CORAL (fork `smith6jt-cop/CORAL` @
  `2aa29e6d`), on the maintainer's instruction of 2026-09-06 that the granted Hugging Face access is for
  KRONOS2 and that it is to be used through the CORAL GitHub. Arm F `kronos-clusters` (KRONOS v1) stays
  `Unavailable(huggingface_access_pending)`. Configuration frozen in the module docstring before its first
  record: CORAL ingest of each frame window (every identity marker on the donor + DAPI), imported frame
  label masks (global labels preserved), Otsu tissue mask, `PatchConfig(64 px, cell_centered)` (32 um at
  0.5078 um/px; CORAL's cell-tutorial tile), cell-isolated tiles, fp32 CLS, batch 128; seven panel
  markers novel to the KRONOS2 vocabulary (Ker8_18, Keratin_5, Caveolin, TP63, GCG, INS, SST) registered
  through CORAL's `additional_markers.csv` with compartment/family reused from existing KRONOS2
  categories (recorded assumptions) and (mean, std) computed by CORAL from the donor's own tissue-masked
  windows; same-antigen renames CD209→DC-SIGN, SMA→A-SMA, M2Gal3→GALECTIN-3, Pan_Cytokeratin→CYTOKERATIN,
  Iba1→IBA1, CD11c→CD11C, CD79a→CD79A, E_cadherin→E-CADHERIN, EpCAM→EPCAM, Vimentin→VIMENTIN,
  Podoplanin→PODOPLANIN. Clusters and annotation identical to arm F (kNN k = 15, Leiden at
  {0.3, 0.5, 1.0} × 5 seeds, null-gated annotation on the marker values, resolution by seed stability);
  the annotation-only false-award rate on a column-shuffled marker block is recorded per run. External
  pins counted for the A.3 tie-break: KRONOS2 weights, KRONOS2 marker metadata, CORAL commit.
