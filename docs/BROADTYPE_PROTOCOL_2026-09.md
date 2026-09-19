---
title: Broad-typing bake-off protocol (PhenoCycler)
status: pre-registered 2026-09-05; Sections 0-3 are frozen before the first screen record
supersedes: none (the DeepCellTypes review of 2026-08-24 stays closed and is cited, not reopened)
creates_cell_calls: development only (data/development/broadtype/); no production output
production_outputs_written: none; registry latch NO_VERIFIED_TYPING unchanged
owner: smith6jt
lane: core/phenocycler/scripts/broadtype/
---

> **STATUS, 2026-09-18 — this comparison is retired and nothing here is live.** The producer choice
> was withdrawn on 2026-09-08 for the twelve verified defects recorded in §4.9 below, and the
> comparison itself is now on the never-revive list in the root `CLAUDE.md` alongside RESTORE, the
> Gate 2-5 apparatus and the GMM/Otsu gate searches. `data/development/broadtype/release/` still
> exists on disk and **must not be read as a verified typing**; `NO_VERIFIED_TYPING` never moved.
> The PhenoCycler typing that stands today is the **subpanel classification** — context from QuPath
> objects, then one question per context against per-donor reference levels drawn from the ordinary
> tissue — in `data/development/typing/cells/`, mechanics in
> `core/phenocycler/scripts/typing/RECONSTRUCTION.md`. Read this document as the frozen record of a
> closed comparison. Nothing above §4's results line is edited, per its own rule; the title's
> "bake-off" and the per-method "arm" identifiers are the vocabulary of the day and are not to be
> carried into new writing — say *method*, or *the comparison*.

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

## Screen summary, appended 2026-09-07 17:26 UTC

Scope: this arm at this configuration on donors 6436 and 6442. Frame: `frame`. Declarations `broadtype-declarations-2026-09-05`; 7 declared labels; survival needs >= 5 verified per donor on both donors (criteria `survive_frac` = 0.6667). Verdict `broadtype-verdict-2026-09-05`, written 2026-09-07T17:26:08Z.

### Per-label verdicts (arm x donor x label)

| arm | donor | label | status | blocking reasons |
|---|---|---|---|---|
| celesta | 6436 | Endocrine | verified | - |
| celesta | 6436 | Epithelial_nonendocrine | failed | heldout: 0/6 markers clear, need 4 (E_cadherin: h=0.530 <= p95=1.485; EpCAM: h=1.426 <= p95=1.711; Ker8_18: h=1.570 <= p95=1.855; Keratin_5: h=1.088 <= p95=1.216; Pan_Cytokeratin: h=0.965 <= p95=1.653; TP63: h=0.902 <= p95=1.424) |
| celesta | 6436 | Vascular | verified | - |
| celesta | 6436 | Stromal | verified | - |
| celesta | 6436 | T_NK | failed | heldout: 2/4 markers clear, need 3 (CD3e: h=1.354 <= p95=2.031; FOXP3: h=1.953 <= p95=1.953) |
| celesta | 6436 | B_Plasma | failed | heldout: 0/2 markers clear, need 2 (CD20: h=1.988 <= p95=1.988; CD79a: h=1.453 <= p95=1.453) |
| celesta | 6436 | Myeloid | verified | - |
| celesta | 6442 | Endocrine | verified | - |
| celesta | 6442 | Epithelial_nonendocrine | failed | heldout: 1/6 markers clear, need 4 (E_cadherin: h=0.676 <= p95=1.201; EpCAM: h=0.711 <= p95=1.189; Ker8_18: h=0.479 <= p95=1.117; Keratin_5: h=1.199 <= p95=1.199; Pan_Cytokeratin: h=0.719 <= p95=1.123) |
| celesta | 6442 | Vascular | verified | - |
| celesta | 6442 | Stromal | failed | heldout: 1/2 markers clear, need 2 (Vimentin: h=0.854 <= p95=1.707) |
| celesta | 6442 | T_NK | verified | - |
| celesta | 6442 | B_Plasma | verified | - |
| celesta | 6442 | Myeloid | failed | heldout: 5/8 markers clear, need 6 (CD206: h=1.653 <= p95=2.522; CD209: h=0.000 <= p95=2.957; CD68: h=0.000 <= p95=3.528) |
| cosine-37-collapsed | 6436 | Endocrine | failed | state_marker_in_identity_inputs: ['Bcl_2', 'Beta_actin', 'CD107a', 'CD38', 'CD39', 'CD44', 'CD57', 'CD99', 'Granzyme_B', 'HLA_A', 'HLA_DR', 'ICOS', 'IDO1', 'Ki67', 'LAG3', 'PCNA', 'PD_1', 'PD_L1', 'SOX2', 'TCF_1', 'TOX', 'VISTA', 'iNOS'] |
| cosine-37-collapsed | 6436 | Epithelial_nonendocrine | failed | state_marker_in_identity_inputs: ['Bcl_2', 'Beta_actin', 'CD107a', 'CD38', 'CD39', 'CD44', 'CD57', 'CD99', 'Granzyme_B', 'HLA_A', 'HLA_DR', 'ICOS', 'IDO1', 'Ki67', 'LAG3', 'PCNA', 'PD_1', 'PD_L1', 'SOX2', 'TCF_1', 'TOX', 'VISTA', 'iNOS'] |
| cosine-37-collapsed | 6436 | Vascular | failed | state_marker_in_identity_inputs: ['Bcl_2', 'Beta_actin', 'CD107a', 'CD38', 'CD39', 'CD44', 'CD57', 'CD99', 'Granzyme_B', 'HLA_A', 'HLA_DR', 'ICOS', 'IDO1', 'Ki67', 'LAG3', 'PCNA', 'PD_1', 'PD_L1', 'SOX2', 'TCF_1', 'TOX', 'VISTA', 'iNOS'] |
| cosine-37-collapsed | 6436 | Stromal | failed | state_marker_in_identity_inputs: ['Bcl_2', 'Beta_actin', 'CD107a', 'CD38', 'CD39', 'CD44', 'CD57', 'CD99', 'Granzyme_B', 'HLA_A', 'HLA_DR', 'ICOS', 'IDO1', 'Ki67', 'LAG3', 'PCNA', 'PD_1', 'PD_L1', 'SOX2', 'TCF_1', 'TOX', 'VISTA', 'iNOS'] |
| cosine-37-collapsed | 6436 | T_NK | failed | state_marker_in_identity_inputs: ['Bcl_2', 'Beta_actin', 'CD107a', 'CD38', 'CD39', 'CD44', 'CD57', 'CD99', 'Granzyme_B', 'HLA_A', 'HLA_DR', 'ICOS', 'IDO1', 'Ki67', 'LAG3', 'PCNA', 'PD_1', 'PD_L1', 'SOX2', 'TCF_1', 'TOX', 'VISTA', 'iNOS'] |
| cosine-37-collapsed | 6436 | B_Plasma | failed | state_marker_in_identity_inputs: ['Bcl_2', 'Beta_actin', 'CD107a', 'CD38', 'CD39', 'CD44', 'CD57', 'CD99', 'Granzyme_B', 'HLA_A', 'HLA_DR', 'ICOS', 'IDO1', 'Ki67', 'LAG3', 'PCNA', 'PD_1', 'PD_L1', 'SOX2', 'TCF_1', 'TOX', 'VISTA', 'iNOS'] |
| cosine-37-collapsed | 6436 | Myeloid | failed | state_marker_in_identity_inputs: ['Bcl_2', 'Beta_actin', 'CD107a', 'CD38', 'CD39', 'CD44', 'CD57', 'CD99', 'Granzyme_B', 'HLA_A', 'HLA_DR', 'ICOS', 'IDO1', 'Ki67', 'LAG3', 'PCNA', 'PD_1', 'PD_L1', 'SOX2', 'TCF_1', 'TOX', 'VISTA', 'iNOS'] |
| cosine-37-collapsed | 6442 | Endocrine | failed | state_marker_in_identity_inputs: ['Bcl_2', 'Beta_actin', 'CD107a', 'CD38', 'CD39', 'CD44', 'CD57', 'CD99', 'Granzyme_B', 'HLA_A', 'HLA_DR', 'ICOS', 'IDO1', 'Ki67', 'LAG3', 'PCNA', 'PD_1', 'PD_L1', 'SOX2', 'TCF_1', 'TOX', 'VISTA', 'iNOS'] |
| cosine-37-collapsed | 6442 | Epithelial_nonendocrine | failed | state_marker_in_identity_inputs: ['Bcl_2', 'Beta_actin', 'CD107a', 'CD38', 'CD39', 'CD44', 'CD57', 'CD99', 'Granzyme_B', 'HLA_A', 'HLA_DR', 'ICOS', 'IDO1', 'Ki67', 'LAG3', 'PCNA', 'PD_1', 'PD_L1', 'SOX2', 'TCF_1', 'TOX', 'VISTA', 'iNOS'] |
| cosine-37-collapsed | 6442 | Vascular | failed | state_marker_in_identity_inputs: ['Bcl_2', 'Beta_actin', 'CD107a', 'CD38', 'CD39', 'CD44', 'CD57', 'CD99', 'Granzyme_B', 'HLA_A', 'HLA_DR', 'ICOS', 'IDO1', 'Ki67', 'LAG3', 'PCNA', 'PD_1', 'PD_L1', 'SOX2', 'TCF_1', 'TOX', 'VISTA', 'iNOS'] |
| cosine-37-collapsed | 6442 | Stromal | failed | state_marker_in_identity_inputs: ['Bcl_2', 'Beta_actin', 'CD107a', 'CD38', 'CD39', 'CD44', 'CD57', 'CD99', 'Granzyme_B', 'HLA_A', 'HLA_DR', 'ICOS', 'IDO1', 'Ki67', 'LAG3', 'PCNA', 'PD_1', 'PD_L1', 'SOX2', 'TCF_1', 'TOX', 'VISTA', 'iNOS'] |
| cosine-37-collapsed | 6442 | T_NK | failed | state_marker_in_identity_inputs: ['Bcl_2', 'Beta_actin', 'CD107a', 'CD38', 'CD39', 'CD44', 'CD57', 'CD99', 'Granzyme_B', 'HLA_A', 'HLA_DR', 'ICOS', 'IDO1', 'Ki67', 'LAG3', 'PCNA', 'PD_1', 'PD_L1', 'SOX2', 'TCF_1', 'TOX', 'VISTA', 'iNOS'] |
| cosine-37-collapsed | 6442 | B_Plasma | failed | state_marker_in_identity_inputs: ['Bcl_2', 'Beta_actin', 'CD107a', 'CD38', 'CD39', 'CD44', 'CD57', 'CD99', 'Granzyme_B', 'HLA_A', 'HLA_DR', 'ICOS', 'IDO1', 'Ki67', 'LAG3', 'PCNA', 'PD_1', 'PD_L1', 'SOX2', 'TCF_1', 'TOX', 'VISTA', 'iNOS'] |
| cosine-37-collapsed | 6442 | Myeloid | failed | state_marker_in_identity_inputs: ['Bcl_2', 'Beta_actin', 'CD107a', 'CD38', 'CD39', 'CD44', 'CD57', 'CD99', 'Granzyme_B', 'HLA_A', 'HLA_DR', 'ICOS', 'IDO1', 'Ki67', 'LAG3', 'PCNA', 'PD_1', 'PD_L1', 'SOX2', 'TCF_1', 'TOX', 'VISTA', 'iNOS'] |
| cosine-broad8 | 6436 | Endocrine | verified | - |
| cosine-broad8 | 6436 | Epithelial_nonendocrine | failed | heldout: 2/6 markers clear, need 4 (E_cadherin: h=1.345 <= p95=1.479; EpCAM: h=0.339 <= p95=1.582; Ker8_18: h=1.039 <= p95=1.336; Keratin_5: h=0.472 <= p95=1.422) |
| cosine-broad8 | 6436 | Vascular | verified | - |
| cosine-broad8 | 6436 | Stromal | failed | heldout: 0/2 markers clear, need 2 (SMA: h=3.598 <= p95=5.397; Vimentin: h=0.000 <= p95=8.858) |
| cosine-broad8 | 6436 | T_NK | verified | - |
| cosine-broad8 | 6436 | B_Plasma | failed | heldout: 1/2 markers clear, need 2 (CD79a: h=1.357 <= p95=2.375) |
| cosine-broad8 | 6436 | Myeloid | verified | - |
| cosine-broad8 | 6442 | Endocrine | verified | - |
| cosine-broad8 | 6442 | Epithelial_nonendocrine | verified | - |
| cosine-broad8 | 6442 | Vascular | verified | - |
| cosine-broad8 | 6442 | Stromal | failed | heldout: 1/2 markers clear, need 2 (Vimentin: h=0.000 <= p95=7.054) |
| cosine-broad8 | 6442 | T_NK | verified | - |
| cosine-broad8 | 6442 | B_Plasma | verified | - |
| cosine-broad8 | 6442 | Myeloid | verified | - |
| dct-ptft-broad8-typeonly | 6436 | Endocrine | unjudgeable | spatial: Endocrine: calls: only 0 usable observations, need >= 50 |
| dct-ptft-broad8-typeonly | 6436 | Epithelial_nonendocrine | verified | - |
| dct-ptft-broad8-typeonly | 6436 | Vascular | failed | heldout: 1/3 markers clear, need 2 (CD31: h=0.000 <= p95=6.465; Podoplanin: h=0.000 <= p95=2.502) |
| dct-ptft-broad8-typeonly | 6436 | Stromal | unjudgeable | heldout: 1 judgeable marker(s), need >= 2 (SMA: Stromal: calls: only 0 usable observations, need >= 50) |
| dct-ptft-broad8-typeonly | 6436 | T_NK | verified | - |
| dct-ptft-broad8-typeonly | 6436 | B_Plasma | failed | heldout: 0/2 markers clear, need 2 (CD20: h=6.502 <= p95=6.502; CD79a: h=7.579 <= p95=7.579) |
| dct-ptft-broad8-typeonly | 6436 | Myeloid | verified | - |
| dct-ptft-broad8-typeonly | 6442 | Endocrine | unjudgeable | spatial: Endocrine: calls: only 0 usable observations, need >= 50 |
| dct-ptft-broad8-typeonly | 6442 | Epithelial_nonendocrine | failed | heldout: 0/6 markers clear, need 4 (E_cadherin: h=1.038 <= p95=1.038; EpCAM: h=1.047 <= p95=1.047; Ker8_18: h=1.011 <= p95=1.036; Keratin_5: h=1.046 <= p95=1.046; Pan_Cytokeratin: h=1.034 <= p95=1.034; TP63: h=0.000 <= p95=2.467) |
| dct-ptft-broad8-typeonly | 6442 | Vascular | failed | heldout: 1/2 markers clear, need 2 (Podoplanin: h=1.746 <= p95=3.579); heldout (not judged): CD31: Vascular: calls: only 23 usable observations, need >= 50 |
| dct-ptft-broad8-typeonly | 6442 | Stromal | unjudgeable | spatial: Stromal: calls: only 20 usable observations, need >= 50; heldout: 0 judgeable marker(s), need >= 2 (SMA: Stromal: calls: only 0 usable observations, need >= 50; Vimentin: Stromal: calls: only 28 usable observations, need >= 50) |
| dct-ptft-broad8-typeonly | 6442 | T_NK | verified | - |
| dct-ptft-broad8-typeonly | 6442 | B_Plasma | unjudgeable | spatial: B_Plasma: calls: only 28 usable observations, need >= 50; heldout: 0 judgeable marker(s), need >= 2 (CD20: B_Plasma: calls: only 1 usable observations, need >= 50; CD79a: B_Plasma: calls: only 16 usable observations, need >= 50) |
| dct-ptft-broad8-typeonly | 6442 | Myeloid | failed | heldout: 4/8 markers clear, need 6 (CD11c: h=5.293 <= p95=5.293; CD163: h=5.765 <= p95=5.765; Iba1: h=6.506 <= p95=6.506; MPO: h=0.000 <= p95=5.911) |
| nimbus-rules | 6436 | Endocrine | verified | - |
| nimbus-rules | 6436 | Epithelial_nonendocrine | verified | - |
| nimbus-rules | 6436 | Vascular | failed | stability: spatial verdict flips on ['subsample:1'] |
| nimbus-rules | 6436 | Stromal | unjudgeable | heldout: 1 judgeable marker(s), need >= 2 (Vimentin: Stromal: calls: only 39 usable observations, need >= 50) |
| nimbus-rules | 6436 | T_NK | verified | - |
| nimbus-rules | 6436 | B_Plasma | failed | spatial: I=0.003 <= p95=0.004 (n_calls=102); heldout: 1 judgeable marker(s), need >= 2 (CD20: B_Plasma: calls: only 42 usable observations, need >= 50) |
| nimbus-rules | 6436 | Myeloid | verified | - |
| nimbus-rules | 6442 | Endocrine | verified | - |
| nimbus-rules | 6442 | Epithelial_nonendocrine | failed | heldout: 3/6 markers clear, need 4 (Ker8_18: h=1.205 <= p95=1.246; Pan_Cytokeratin: h=0.976 <= p95=1.249; TP63: h=1.067 <= p95=1.186) |
| nimbus-rules | 6442 | Vascular | failed | heldout: 0/3 markers clear, need 2 (CD31: h=0.000 <= p95=0.000; Caveolin: h=0.000 <= p95=0.633; Podoplanin: h=0.000 <= p95=11.306) |
| nimbus-rules | 6442 | Stromal | failed | heldout: 0/2 markers clear, need 2 (SMA: h=6.692 <= p95=6.692; Vimentin: h=0.000 <= p95=1.040) |
| nimbus-rules | 6442 | T_NK | verified | - |
| nimbus-rules | 6442 | B_Plasma | failed | stability: spatial verdict flips on ['subsample:3', 'subsample:4']; heldout: 1 judgeable marker(s), need >= 2 (CD20: B_Plasma: calls: only 2 usable observations, need >= 50) |
| nimbus-rules | 6442 | Myeloid | verified | - |
| som-kintsugi | 6436 | Endocrine | unjudgeable | spatial: Endocrine: calls: only 0 usable observations, need >= 50 |
| som-kintsugi | 6436 | Epithelial_nonendocrine | failed | stability: spatial verdict flips on ['seed:1']; spatial: Epithelial_nonendocrine: calls: only 0 usable observations, need >= 50; heldout: 0 judgeable marker(s), need >= 2 (E_cadherin: Epithelial_nonendocrine: tail cells in pool: only 0 usable observations, need >= 20; EpCAM: Epithelial_nonendocrine: tail cells in pool: only 0 usable observations, need >= 20; Ker8_18: Epithelial_nonendocrine: tail cells in pool: only 0 usable observations, need >= 20; Keratin_5: Epithelial_nonendocrine: tail cells in pool: only 16 usable observations, need >= 20; Pan_Cytokeratin: Epithelial_nonendocrine: tail cells in pool: only 0 usable observations, need >= 20; TP63: Epithelial_nonendocrine: tail cells in pool: only 19 usable observations, need >= 20) |
| som-kintsugi | 6436 | Vascular | unjudgeable | spatial: Vascular: calls: only 0 usable observations, need >= 50; heldout: 0 judgeable marker(s), need >= 2 (CD31: Vascular: calls: only 0 usable observations, need >= 50; Caveolin: Vascular: calls: only 0 usable observations, need >= 50; Podoplanin: Vascular: calls: only 0 usable observations, need >= 50) |
| som-kintsugi | 6436 | Stromal | failed | stability: spatial verdict flips on ['seed:1']; spatial: Stromal: calls: only 0 usable observations, need >= 50; heldout: 0 judgeable marker(s), need >= 2 (SMA: Stromal: calls: only 0 usable observations, need >= 50; Vimentin: Stromal: calls: only 0 usable observations, need >= 50) |
| som-kintsugi | 6436 | T_NK | failed | stability: spatial verdict flips on ['seed:1', 'seed:3']; spatial: T_NK: calls: only 0 usable observations, need >= 50; heldout: 1 judgeable marker(s), need >= 2 (CD4: T_NK: calls: only 0 usable observations, need >= 50; CD8: T_NK: tail cells in pool: only 15 usable observations, need >= 20; FOXP3: T_NK: calls: only 0 usable observations, need >= 50) |
| som-kintsugi | 6436 | B_Plasma | unjudgeable | spatial: B_Plasma: calls: only 0 usable observations, need >= 50; heldout: 0 judgeable marker(s), need >= 2 (CD20: B_Plasma: calls: only 0 usable observations, need >= 50; CD79a: B_Plasma: calls: only 0 usable observations, need >= 50) |
| som-kintsugi | 6436 | Myeloid | unjudgeable | spatial: Myeloid: calls: only 0 usable observations, need >= 50; heldout: 0 judgeable marker(s), need >= 2 (CD11c: Myeloid: calls: only 0 usable observations, need >= 50; CD163: Myeloid: calls: only 0 usable observations, need >= 50; CD206: Myeloid: calls: only 0 usable observations, need >= 50; CD209: Myeloid: calls: only 0 usable observations, need >= 50; CD68: Myeloid: calls: only 0 usable observations, need >= 50; Iba1: Myeloid: tail cells in pool: only 10 usable observations, need >= 20; M2Gal3: Myeloid: tail cells in pool: only 15 usable observations, need >= 20; MPO: Myeloid: calls: only 0 usable observations, need >= 50) |
| som-kintsugi | 6442 | Endocrine | failed | stability: spatial verdict flips on ['seed:3', 'seed:4']; spatial: Endocrine: calls: only 0 usable observations, need >= 50 |
| som-kintsugi | 6442 | Epithelial_nonendocrine | failed | heldout: 1/4 markers clear, need 3 (EpCAM: h=0.000 <= p95=1.375; Keratin_5: h=0.765 <= p95=1.529; TP63: h=0.000 <= p95=1.167); stability: spatial verdict flips on ['seed:3', 'seed:4']; spatial: Epithelial_nonendocrine: calls: only 0 usable observations, need >= 50; heldout (not judged): E_cadherin: Epithelial_nonendocrine: tail cells in pool: only 0 usable observations, need >= 20; Pan_Cytokeratin: Epithelial_nonendocrine: calls: only 0 usable observations, need >= 50 |
| som-kintsugi | 6442 | Vascular | unjudgeable | spatial: Vascular: calls: only 0 usable observations, need >= 50; heldout: 0 judgeable marker(s), need >= 2 (CD31: Vascular: calls: only 0 usable observations, need >= 50; Caveolin: Vascular: tail cells in pool: only 18 usable observations, need >= 20; Podoplanin: Vascular: calls: only 0 usable observations, need >= 50) |
| som-kintsugi | 6442 | Stromal | unjudgeable | heldout: 0 judgeable marker(s), need >= 2 (SMA: Stromal: calls: only 0 usable observations, need >= 50; Vimentin: Stromal: calls: only 0 usable observations, need >= 50) |
| som-kintsugi | 6442 | T_NK | unjudgeable | spatial: T_NK: calls: only 0 usable observations, need >= 50; heldout: 0 judgeable marker(s), need >= 2 (CD3e: T_NK: tail cells in pool: only 19 usable observations, need >= 20; CD4: T_NK: tail cells in pool: only 8 usable observations, need >= 20; CD8: T_NK: calls: only 0 usable observations, need >= 50; FOXP3: T_NK: calls: only 0 usable observations, need >= 50) |
| som-kintsugi | 6442 | B_Plasma | unjudgeable | spatial: B_Plasma: calls: only 0 usable observations, need >= 50; heldout: 0 judgeable marker(s), need >= 2 (CD20: B_Plasma: tail cells in pool: only 19 usable observations, need >= 20; CD79a: B_Plasma: calls: only 0 usable observations, need >= 50) |
| som-kintsugi | 6442 | Myeloid | unjudgeable | spatial: Myeloid: calls: only 0 usable observations, need >= 50; heldout: 0 judgeable marker(s), need >= 2 (CD11c: Myeloid: calls: only 0 usable observations, need >= 50; CD163: Myeloid: calls: only 0 usable observations, need >= 50; CD206: Myeloid: tail cells in pool: only 7 usable observations, need >= 20; CD209: Myeloid: calls: only 0 usable observations, need >= 50; CD68: Myeloid: calls: only 0 usable observations, need >= 50; Iba1: Myeloid: calls: only 0 usable observations, need >= 50; M2Gal3: Myeloid: calls: only 0 usable observations, need >= 50; MPO: Myeloid: calls: only 0 usable observations, need >= 50) |
| som-kintsugi-aligned | 6436 | Endocrine | verified | - |
| som-kintsugi-aligned | 6436 | Epithelial_nonendocrine | failed | stability: spatial verdict flips on ['seed:1', 'seed:3']; spatial: Epithelial_nonendocrine: calls: only 0 usable observations, need >= 50; heldout: 1 judgeable marker(s), need >= 2 (E_cadherin: Epithelial_nonendocrine: calls: only 0 usable observations, need >= 50; EpCAM: Epithelial_nonendocrine: calls: only 0 usable observations, need >= 50; Ker8_18: Epithelial_nonendocrine: calls: only 0 usable observations, need >= 50; Pan_Cytokeratin: Epithelial_nonendocrine: calls: only 0 usable observations, need >= 50; TP63: Epithelial_nonendocrine: tail cells in pool: only 17 usable observations, need >= 20) |
| som-kintsugi-aligned | 6436 | Vascular | failed | stability: spatial verdict flips on ['seed:4']; heldout: 0 judgeable marker(s), need >= 2 (CD31: Vascular: tail cells in pool: only 5 usable observations, need >= 20; Caveolin: Vascular: tail cells in pool: only 3 usable observations, need >= 20; Podoplanin: Vascular: tail cells in pool: only 16 usable observations, need >= 20) |
| som-kintsugi-aligned | 6436 | Stromal | unjudgeable | heldout: 0 judgeable marker(s), need >= 2 (SMA: Stromal: tail cells in pool: only 13 usable observations, need >= 20; Vimentin: Stromal: calls: only 0 usable observations, need >= 50) |
| som-kintsugi-aligned | 6436 | T_NK | verified | - |
| som-kintsugi-aligned | 6436 | B_Plasma | unjudgeable | spatial: B_Plasma: calls: only 0 usable observations, need >= 50; heldout: 0 judgeable marker(s), need >= 2 (CD20: B_Plasma: tail cells in pool: only 2 usable observations, need >= 20; CD79a: B_Plasma: calls: only 0 usable observations, need >= 50) |
| som-kintsugi-aligned | 6436 | Myeloid | failed | heldout: 3/6 markers clear, need 4 (CD68: h=5.260 <= p95=5.260; M2Gal3: h=0.000 <= p95=1.581; MPO: h=14.048 <= p95=14.048); heldout (not judged): CD11c: Myeloid: calls: only 0 usable observations, need >= 50; CD163: Myeloid: calls: only 0 usable observations, need >= 50 |
| som-kintsugi-aligned | 6442 | Endocrine | verified | - |
| som-kintsugi-aligned | 6442 | Epithelial_nonendocrine | failed | heldout: 0/2 markers clear, need 2 (EpCAM: h=1.009 <= p95=1.345; Keratin_5: h=0.849 <= p95=1.273); stability: spatial verdict flips on ['seed:3']; heldout (not judged): E_cadherin: Epithelial_nonendocrine: calls: only 0 usable observations, need >= 50; Ker8_18: Epithelial_nonendocrine: calls: only 0 usable observations, need >= 50; Pan_Cytokeratin: Epithelial_nonendocrine: calls: only 0 usable observations, need >= 50; TP63: Epithelial_nonendocrine: tail cells in pool: only 8 usable observations, need >= 20 |
| som-kintsugi-aligned | 6442 | Vascular | failed | stability: spatial verdict flips on ['seed:1']; spatial: Vascular: calls: only 0 usable observations, need >= 50; heldout: 1 judgeable marker(s), need >= 2 (CD31: Vascular: tail cells in pool: only 7 usable observations, need >= 20; Caveolin: Vascular: tail cells in pool: only 3 usable observations, need >= 20) |
| som-kintsugi-aligned | 6442 | Stromal | unjudgeable | heldout: 0 judgeable marker(s), need >= 2 (SMA: Stromal: calls: only 0 usable observations, need >= 50; Vimentin: Stromal: calls: only 0 usable observations, need >= 50) |
| som-kintsugi-aligned | 6442 | T_NK | verified | - |
| som-kintsugi-aligned | 6442 | B_Plasma | unjudgeable | spatial: B_Plasma: calls: only 0 usable observations, need >= 50; heldout: 0 judgeable marker(s), need >= 2 (CD20: B_Plasma: tail cells in pool: only 15 usable observations, need >= 20; CD79a: B_Plasma: calls: only 0 usable observations, need >= 50) |
| som-kintsugi-aligned | 6442 | Myeloid | failed | heldout: 4/7 markers clear, need 5 (CD11c: h=0.000 <= p95=5.286; Iba1: h=0.000 <= p95=2.524; MPO: h=0.000 <= p95=4.965); stability: spatial verdict flips on ['seed:3', 'seed:4', 'seed:5']; heldout (not judged): M2Gal3: Myeloid: calls: only 0 usable observations, need >= 50 |
| tacit | 6436 | Endocrine | verified | - |
| tacit | 6436 | Epithelial_nonendocrine | verified | - |
| tacit | 6436 | Vascular | verified | - |
| tacit | 6436 | Stromal | failed | heldout: 0/2 markers clear, need 2 (SMA: h=0.856 <= p95=1.712; Vimentin: h=0.000 <= p95=2.285) |
| tacit | 6436 | T_NK | failed | heldout: 1/4 markers clear, need 3 (CD3e: h=0.459 <= p95=1.838; CD4: h=0.864 <= p95=2.593; CD8: h=0.484 <= p95=1.936); stability: spatial verdict flips on ['subsample:1'] |
| tacit | 6436 | B_Plasma | failed | spatial: I=0.012 <= p95=0.012 (n_calls=636); heldout: 0/2 markers clear, need 2 (CD20: h=0.000 <= p95=4.788; CD79a: h=0.577 <= p95=1.732); stability: spatial verdict flips on ['subsample:1', 'subsample:5'] |
| tacit | 6436 | Myeloid | verified | - |
| tacit | 6442 | Endocrine | verified | - |
| tacit | 6442 | Epithelial_nonendocrine | verified | - |
| tacit | 6442 | Vascular | failed | heldout: 0/2 markers clear, need 2 (CD31: h=0.000 <= p95=16.907; Podoplanin: h=0.000 <= p95=7.685); stability: spatial verdict flips on ['subsample:1', 'subsample:3', 'subsample:4']; heldout (not judged): Caveolin: Vascular: calls: only 28 usable observations, need >= 50 |
| tacit | 6442 | Stromal | failed | heldout: 1/2 markers clear, need 2 (Vimentin: h=0.701 <= p95=2.138) |
| tacit | 6442 | T_NK | failed | heldout: 0/2 markers clear, need 2 (CD3e: h=0.000 <= p95=2.082; CD4: h=0.000 <= p95=5.149); stability: spatial verdict flips on ['subsample:1', 'subsample:4', 'subsample:5']; heldout (not judged): CD8: T_NK: calls: only 20 usable observations, need >= 50; FOXP3: T_NK: calls: only 34 usable observations, need >= 50 |
| tacit | 6442 | B_Plasma | failed | heldout: 1/2 markers clear, need 2 (CD20: h=0.321 <= p95=1.285); stability: spatial verdict flips on ['subsample:1', 'subsample:3', 'subsample:4', 'subsample:5'] |
| tacit | 6442 | Myeloid | failed | heldout: 2/8 markers clear, need 6 (CD11c: h=9.848 <= p95=9.848; CD163: h=4.090 <= p95=4.090; CD68: h=4.718 <= p95=4.718; Iba1: h=0.000 <= p95=13.896; M2Gal3: h=0.000 <= p95=3.857; MPO: h=0.000 <= p95=10.567) |
| tasp | 6436 | Endocrine | verified | - |
| tasp | 6436 | Epithelial_nonendocrine | verified | - |
| tasp | 6436 | Vascular | verified | - |
| tasp | 6436 | Stromal | failed | heldout: 0/2 markers clear, need 2 (SMA: h=0.000 <= p95=17.378; Vimentin: h=0.000 <= p95=10.233) |
| tasp | 6436 | T_NK | verified | - |
| tasp | 6436 | B_Plasma | unjudgeable | heldout: 1 judgeable marker(s), need >= 2 (CD20: B_Plasma: calls: only 20 usable observations, need >= 50) |
| tasp | 6436 | Myeloid | verified | - |
| tasp | 6442 | Endocrine | verified | - |
| tasp | 6442 | Epithelial_nonendocrine | verified | - |
| tasp | 6442 | Vascular | verified | - |
| tasp | 6442 | Stromal | failed | heldout: 0/2 markers clear, need 2 (SMA: h=0.552 <= p95=1.545; Vimentin: h=0.000 <= p95=4.280) |
| tasp | 6442 | T_NK | verified | - |
| tasp | 6442 | B_Plasma | unjudgeable | spatial: B_Plasma: calls: only 42 usable observations, need >= 50; heldout: 0 judgeable marker(s), need >= 2 (CD20: B_Plasma: calls: only 31 usable observations, need >= 50; CD79a: B_Plasma: calls: only 41 usable observations, need >= 50) |
| tasp | 6442 | Myeloid | verified | - |

Scope of every row above: this arm at this configuration on donors 6436 and 6442. Noise false-award rate measured on synthetic noise (from the arms' records): som-kintsugi: 6436: {"false_award_measurement": {"shuffled_copy": "every feature column independently permuted (rng seed ctx.seed), then standardised", "annotation_only": "annotation on the shuffled Z with the PRIMARY clusters; awarded / clusters >= min_cells (the gate's Type-I rate on random cell sets)", "arm": "SOM re-fitted on the shuffled block over the same nodes x seeds grid, annotated with the same gate; awarded / clusters >= min_cells pooled over the 5 seeds per node count (what the arm does on structureless input)", "headline": "false_award_rate = pooled arm rate at the chosen node count"}}; 6442: {"false_award_measurement": {"shuffled_copy": "every feature column independently permuted (rng seed ctx.seed), then standardised", "annotation_only": "annotation on the shuffled Z with the PRIMARY clusters; awarded / clusters >= min_cells (the gate's Type-I rate on random cell sets)", "arm": "SOM re-fitted on the shuffled block over the same nodes x seeds grid, annotated with the same gate; awarded / clusters >= min_cells pooled over the 5 seeds per node count (what the arm does on structureless input)", "headline": "false_award_rate = pooled arm rate at the chosen node count"}} | som-kintsugi-aligned: 6436: {"false_award_measurement": {"shuffled_copy": "every feature column independently permuted (rng seed ctx.seed), then standardised", "annotation_only": "annotation on the shuffled Z with the PRIMARY clusters; awarded / clusters >= min_cells (the gate's Type-I rate on random cell sets)", "arm": "SOM re-fitted on the shuffled block over the same nodes x seeds grid, annotated with the same gate; awarded / clusters >= min_cells pooled over the 5 seeds per node count (what the arm does on structureless input)", "headline": "false_award_rate = pooled arm rate at the chosen node count"}}; 6442: {"false_award_measurement": {"shuffled_copy": "every feature column independently permuted (rng seed ctx.seed), then standardised", "annotation_only": "annotation on the shuffled Z with the PRIMARY clusters; awarded / clusters >= min_cells (the gate's Type-I rate on random cell sets)", "arm": "SOM re-fitted on the shuffled block over the same nodes x seeds grid, annotated with the same gate; awarded / clusters >= min_cells pooled over the 5 seeds per node count (what the arm does on structureless input)", "headline": "false_award_rate = pooled arm rate at the chosen node count"}}.

### Survival

| arm | 6436 verified | 6442 verified | needed | survives | shortfall | mean held-out margin | external pins | primary wall s |
|---|---|---|---|---|---|---|---|---|
| celesta | 4 | 4 | 5 | no | 6436: 4 verified < 5; 6442: 4 verified < 5 | 4.971 | 0 | 62.8 |
| cosine-37-collapsed | 0 | 0 | 5 | no | 6436: 0 verified < 5; 6442: 0 verified < 5 | nan | 0 | 4.8 |
| cosine-broad8 | 4 | 6 | 5 | no | 6436: 4 verified < 5 | 6.143 | 0 | 12.2 |
| dct-ptft-broad8-typeonly | 3 | 1 | 5 | no | 6436: 3 verified < 5; 6442: 1 verified < 5 | 17.297 | 0 | 524.7 |
| nimbus-rules | 4 | 3 | 5 | no | 6436: 4 verified < 5; 6442: 3 verified < 5 | 13.252 | 0 | 1852.6 |
| som-kintsugi | 0 | 0 | 5 | no | 6436: 0 verified < 5; 6442: 0 verified < 5 | 3.368 | 0 | 166.8 |
| som-kintsugi-aligned | 2 | 2 | 5 | no | 6436: 2 verified < 5; 6442: 2 verified < 5 | 32.025 | 0 | 154.6 |
| tacit | 4 | 2 | 5 | no | 6436: 4 verified < 5; 6442: 2 verified < 5 | 6.782 | 0 | 252.2 |
| tasp | 5 | 5 | 5 | yes | - | 5.382 | 0 | 34.7 |

Scope of the survival row: this arm at this configuration on donors 6436 and 6442. Noise false-award rate measured on synthetic noise (from the arms' records): som-kintsugi: 6436: {"false_award_measurement": {"shuffled_copy": "every feature column independently permuted (rng seed ctx.seed), then standardised", "annotation_only": "annotation on the shuffled Z with the PRIMARY clusters; awarded / clusters >= min_cells (the gate's Type-I rate on random cell sets)", "arm": "SOM re-fitted on the shuffled block over the same nodes x seeds grid, annotated with the same gate; awarded / clusters >= min_cells pooled over the 5 seeds per node count (what the arm does on structureless input)", "headline": "false_award_rate = pooled arm rate at the chosen node count"}}; 6442: {"false_award_measurement": {"shuffled_copy": "every feature column independently permuted (rng seed ctx.seed), then standardised", "annotation_only": "annotation on the shuffled Z with the PRIMARY clusters; awarded / clusters >= min_cells (the gate's Type-I rate on random cell sets)", "arm": "SOM re-fitted on the shuffled block over the same nodes x seeds grid, annotated with the same gate; awarded / clusters >= min_cells pooled over the 5 seeds per node count (what the arm does on structureless input)", "headline": "false_award_rate = pooled arm rate at the chosen node count"}} | som-kintsugi-aligned: 6436: {"false_award_measurement": {"shuffled_copy": "every feature column independently permuted (rng seed ctx.seed), then standardised", "annotation_only": "annotation on the shuffled Z with the PRIMARY clusters; awarded / clusters >= min_cells (the gate's Type-I rate on random cell sets)", "arm": "SOM re-fitted on the shuffled block over the same nodes x seeds grid, annotated with the same gate; awarded / clusters >= min_cells pooled over the 5 seeds per node count (what the arm does on structureless input)", "headline": "false_award_rate = pooled arm rate at the chosen node count"}}; 6442: {"false_award_measurement": {"shuffled_copy": "every feature column independently permuted (rng seed ctx.seed), then standardised", "annotation_only": "annotation on the shuffled Z with the PRIMARY clusters; awarded / clusters >= min_cells (the gate's Type-I rate on random cell sets)", "arm": "SOM re-fitted on the shuffled block over the same nodes x seeds grid, annotated with the same gate; awarded / clusters >= min_cells pooled over the 5 seeds per node count (what the arm does on structureless input)", "headline": "false_award_rate = pooled arm rate at the chosen node count"}}.

### Producer (A.3)

**producer: tasp** -- ranked survivors: tasp [-10, -5.381625675806636, 0, 34.74511456489563]. Rule: survive iff verified >= ceil(survive_frac x n_declared) on BOTH donors; producer = first survivor by (i) verified total, (ii) mean held-out margin over cleared markers, (iii) fewest external pins then lowest primary wall time; agreement never chooses. Scope: this arm at this configuration on donors 6436 and 6442.

## Screen summary, appended 2026-09-07 17:26 UTC

Scope: this arm at this configuration on donors 6436 and 6442. Frame: `donor`. Declarations `broadtype-declarations-2026-09-05`; 7 declared labels; survival needs >= 5 verified per donor on both donors (criteria `survive_frac` = 0.6667). Verdict `broadtype-verdict-2026-09-05`, written 2026-09-07T17:26:09Z.

### Per-label verdicts (arm x donor x label)

| arm | donor | label | status | blocking reasons |
|---|---|---|---|---|
| cosine-broad8 | 6436 | Endocrine | verified | - |
| cosine-broad8 | 6436 | Epithelial_nonendocrine | failed | heldout: 2/6 markers clear, need 4 (E_cadherin: h=0.828 <= p95=1.069; EpCAM: h=0.165 <= p95=1.072; Ker8_18: h=0.621 <= p95=1.075; Pan_Cytokeratin: h=0.933 <= p95=1.079) |
| cosine-broad8 | 6436 | Vascular | verified | - |
| cosine-broad8 | 6436 | Stromal | failed | heldout: 1/2 markers clear, need 2 (Vimentin: h=1.379 <= p95=1.896) |
| cosine-broad8 | 6436 | T_NK | verified | - |
| cosine-broad8 | 6436 | B_Plasma | verified | - |
| cosine-broad8 | 6436 | Myeloid | verified | - |
| cosine-broad8 | 6442 | Endocrine | verified | - |
| cosine-broad8 | 6442 | Epithelial_nonendocrine | verified | - |
| cosine-broad8 | 6442 | Vascular | verified | - |
| cosine-broad8 | 6442 | Stromal | failed | heldout: 1/2 markers clear, need 2 (Vimentin: h=0.444 <= p95=1.778) |
| cosine-broad8 | 6442 | T_NK | verified | - |
| cosine-broad8 | 6442 | B_Plasma | verified | - |
| cosine-broad8 | 6442 | Myeloid | verified | - |
| tasp | 6436 | Endocrine | verified | - |
| tasp | 6436 | Epithelial_nonendocrine | verified | - |
| tasp | 6436 | Vascular | verified | - |
| tasp | 6436 | Stromal | verified | - |
| tasp | 6436 | T_NK | verified | - |
| tasp | 6436 | B_Plasma | verified | - |
| tasp | 6436 | Myeloid | verified | - |
| tasp | 6442 | Endocrine | verified | - |
| tasp | 6442 | Epithelial_nonendocrine | verified | - |
| tasp | 6442 | Vascular | verified | - |
| tasp | 6442 | Stromal | failed | heldout: 1/2 markers clear, need 2 (Vimentin: h=0.817 <= p95=1.838) |
| tasp | 6442 | T_NK | verified | - |
| tasp | 6442 | B_Plasma | verified | - |
| tasp | 6442 | Myeloid | verified | - |

Scope of every row above: this arm at this configuration on donors 6436 and 6442. Noise false-award rate: not recorded in any arm's records (only null-gated cluster arms measure one).

### Survival

| arm | 6436 verified | 6442 verified | needed | survives | shortfall | mean held-out margin | external pins | primary wall s |
|---|---|---|---|---|---|---|---|---|
| cosine-broad8 | 5 | 6 | 5 | yes | - | 6.358 | 0 | 3.4 |
| tasp | 7 | 6 | 5 | yes | - | 13.758 | 0 | 218.2 |

Scope of the survival row: this arm at this configuration on donors 6436 and 6442. Noise false-award rate: not recorded in any arm's records (only null-gated cluster arms measure one).

### Producer (A.3)

**producer: tasp** -- ranked survivors: tasp [-13, -13.758361199942314, 0, 218.2061746120453], cosine-broad8 [-11, -6.358107465660469, 0, 3.4406323432922363]. Rule: survive iff verified >= ceil(survive_frac x n_declared) on BOTH donors; producer = first survivor by (i) verified total, (ii) mean held-out margin over cleared markers, (iii) fewest external pins then lowest primary wall time; agreement never chooses. Scope: this arm at this configuration on donors 6436 and 6442.

## 4.9 Defects found 2026-09-08 — the producer choice is WITHDRAWN

Prompted by the maintainer ("before such a test is reliable, we need to be sure each cell typing method is
being done properly and the data preprocessing is sound"), the lane was audited against the pipeline that
feeds it, each package's own documentation, and the gate's statistics. The audit found defects that
invalidate the producer choice. **No producer is chosen. `producer: tasp` in the tables above is withdrawn.**
The registry latch `NO_VERIFIED_TYPING` was already down and nothing was released on this basis; the
development release `broadtype-tasp-dev-20260906` remains on disk and must not be read as a verified typing.
Every screen record stands as written; nothing above this line is edited.

### Verified by re-reading the code and recomputing from the data

1. **The bake-off compared pipelines, not methods.** Seven arms read the REDSEA + UniFORM parquet; four
   (`dct-ptft-broad8-typeonly`, `nimbus-rules`, `kronos-clusters`, `kronos2-coral`) read raw qptiff pixels in
   the frame windows and never touch REDSEA, UniFORM or the QC flags. Recorded `n_unavailable` on the 6436
   frame, before any method has said anything: 0% (dct-ptft, nimbus), 11.9% (tasp/tacit/celesta), 17.4%
   (cosine-broad8), 72.6% (som-kintsugi), 100% (kronos-clusters).
2. **The cell set is dirty and only one arm cleans it.** `arm_base.load_context` never filters
   `qc_analysis_eligible`; only `cosine_broad8.py:231` does. 17.4% / 16.7% of frame cells are
   analysis-ineligible. 12.9% / 13.1% are NaN in every marker — the objects `redsea.py:674` erases from the
   mask, so their signal is also never subtracted from the neighbours they touch, and the loss is spatially
   structured (9.8%–25.8% across the 16 frame regions of 6436). This lane's own KINTSUGI predecessor filters
   these cells deliberately.
3. **The held-out ablation inverts the evidence it claims to collect.** Withholding a marker a label genuinely
   depends on collapses that label below `min_calls`, which scores as *unjudgeable*, which is not verified.
   Measured: dct-ptft on 6442 drops to 23 Vascular calls without CD31, 1 B_Plasma call without CD20, 0
   Stromal calls without SMA. The labels that verify are the ones whose markers are **redundant**.
4. **`tasp`'s 6436 verdict mixes two versions of the arm.** The `primary` the gate scored for spatial
   coherence and stability was written 2026-09-05T14:01:40Z with `within_donor_rank` features and
   `max_iterations` 20; all 25 ablations and 5 replicates were written 2026-09-06T02:26–02:35Z with
   `standardise_within_donor`, `max_iterations` 500 and a new isolated-seed filter. The held-out gate compared
   new-configuration ablations against a stale-configuration base. Section 1.3's own rule — a configuration
   change after a result is a NEW arm name — was applied to `som-kintsugi-aligned` and not to `tasp`.
5. **`tasp`'s mechanism is the gate's statistic.** It seeds a cell as label L iff the cell is in
   `gate.tail_mask` of a declared marker of L; the gate asks whether L is enriched in `gate.tail_mask` of the
   withheld declared marker. Same function, same q. On 6436, 47–79% of a withheld marker's tail already sits
   inside a co-declared marker's tail for exactly the labels tasp verified.
6. **The tail is not a sane construct on these distributions.** k = 41–42 cells; QC-failed fragments are
   enriched in it up to 7.2x (45% of the 6436 Pan_Cytokeratin tail); tail cells are physically smaller than
   typical, the population REDSEA's perimeter/area term inflates. On the T1D donor the sorted top INS values
   are 7345, 2178, 1874, 1629, then 230, 157 … — four bright cells, then background, so a 42-cell INS tail is
   ~90% background.
7. **Two of the 28 identity markers are a different measurement.** FOXP3 and TP63 are raw uncompensated QuPath
   *nucleus* means; the other 26 are REDSEA-corrected *cell* means. `scripts/celltype/panel.py:63-68` says
   these "must not be compared to the others as though it were the same one". Their dynamic range is ~2x
   median versus ~150–800x for SMA, CD31, GCG.
8. **UniFORM failed its own bar on six identity markers and the failed scales were applied.** `BAR_MODE_SD`
   0.10; Vimentin's residual is 0.852 (47.9x fold across donors) and 6442's Vimentin was still multiplied by
   3.808 against 6436's 1.331. Stromal, whose two declared markers are SMA and Vimentin and which therefore
   needs BOTH held-out tests to clear, failed for essentially every arm.
9. **No multiplicity correction anywhere in the blocking path** (853 null-based tests), and 7.4% of judgeable
   held-out tests were resolved by an exact tie that the strict `>` scores as a failure.
10. **The verdict is not invariant to a cosmetic choice.** `Ambiguous` stays in the statistic's pool,
    `Unavailable` is removed, and nothing constrains which token an arm writes for the same epistemic state.
11. **Arms run in configurations their own authors or this repo would reject.** `dct-ptft-broad8-typeonly` was
    given the `A1-type-only` channel shape that `docs/DEEPCELL_TYPES_REVIEW_2026-08-24.md` had already
    measured on this checkpoint and rejected in bold ("Reducing the panel to bare lineage markers is the wrong
    move"). `som-kintsugi` omits FlowSOM/KINTSUGI's meta-clustering entirely (`som_adapter.meta_cluster`
    exists and is never called). `nimbus-rules` was scored on predictions calling FOXP3 positive in 68% of
    cells and INS in 0.015%. `tacit` is fed unnormalised intensities into a scorer with an absolute internal
    constant. `kronos2-coral` has 6 of ~31 records, no donor 6442, and is absent from `verdict/summary.json`
    — the producer was declared while it was still writing.
12. **The release's provenance pin is broken.** It names the frame summary while being a donor-scope release,
    and the recorded `verdict_summary_sha256` no longer matches that file.

### What still stands

- "This arm produced no usable calls" is a fact where it was recorded (`cosine-37-collapsed`, `som-kintsugi`,
  `kronos-clusters`).
- The mechanical bookkeeping is re-derivable: `verdict.py` recomputes survival from the recorded statistics,
  `gates.require_evidence` fails closed, and the tiering rule is enforced (a record whose identity inputs
  contain a state marker fails every label).
- No composition, disease ordering or expected frequency entered either gate; `assert_frequency_free` is real.
- The QC figures (P01–P05, X01–X03) are measurements of the run and are unaffected by any of the above.
