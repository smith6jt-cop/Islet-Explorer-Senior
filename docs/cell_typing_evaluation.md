# Cell-Typing Workflow Evaluation — Islet Explorer Senior (59-plex CODEX)

_Evaluation + implementation of an automated cell-typing workflow to replace manual `napariGater`
gating. Grounded in a 13-agent research/design/verification pass and empirical tests on donor 6374
(200,000 cells). June 2026._

> **⚠️ STATUS UPDATE (late June 2026) — §§3–8 below are SUPERSEDED.** The global-clustering /
> per-marker-gating / Astir approach evaluated below was abandoned: it produced 25–34% **"Unassigned,"**
> which the user (correctly) rejected as a symptom of **upstream artifacts, not a real category**. The
> current pipeline fixes those artifacts in order **REDSEA → RESTORE → broad-lineage assignment**, then
> per-lineage subclustering. It types **every cell with ZERO "Unassigned."** Full detail in
> `scripts/CLAUDE.md`; summary:
>
> 1. **REDSEA** (`scripts/senior/redsea_full.py`) — scalable pixel-level lateral-spillover correction
>    (the "pan-positive bright" cells were boundary spillover). Validated on 15 donors; impossible
>    cross-compartment co-positives −89%, biology preserved.
> 2. **RESTORE** (`scripts/senior/restore_normalize.py`) — per-image + autofluorescence normalization via
>    mutually-exclusive marker pairs. **Pan_Cytokeratin is the universal negative control** (abundant,
>    immune-negative, bright); the immune **false-positive** bug (bright-acinar autofluorescence tipping a
>    flat threshold) was fixed by **CD20/CD163 ← Pan_Cytokeratin**, with **CD3e ← CD163** kept (Pan_CK
>    destroys CD3e's T-gain). T1D biology preserved & sharpened.
> 3. **Broad lineage** (`scripts/senior/assign_broad_lineage.py`) — hierarchical `_pos` + structural
>    argmax → Epithelial / Fibroblast / Immune / Endocrine / Endothelial / Muscle, 0 Unassigned.
> 4. **Next:** per-lineage Harmony+Leiden subclustering (Step 2); inspect in QuPath via
>    `export_broad_class_for_qupath.py` + `scripts/groovy/import_broad_lineage.groovy` (UUID match).
>
> The §§3–8 evaluation below remains a useful record of *why* the clustering/gating/Astir routes were
> ruled out on this panel (marker quality, batch effects), but is **not** the current method.

## 1. Problem

Phenotyping 23.28 M cells × 59 protein markers (15 donors) currently requires manually dragging
**59 napariGater sliders and re-validating them per donor** before `sm.pp.rescale` →
`sm.tl.phenotype_cells`. This is too tedious to scale and blocks development. The goal: an automated,
multivariate, batch-robust, confidence-scored workflow that confines human effort to **cluster-level
QC (tens of decisions)**, not 59 sliders × 15 donors.

## 2. What the evaluation ruled out (verified on this data)

A research → competing-design → adversarial-verification workflow refuted the obvious "just automate
it" ideas, each with empirical evidence:

| Idea | Verdict |
|---|---|
| Classify on the existing scVI latent | scVI latent encodes **batch, not cell type** (silhouette 0.08 biology vs 0.25 donor). Keep scVI for the trajectory only. |
| Per-donor z-score to remove batch | Corrects mean+scale only; per-donor skew (8–19) / kurtosis (97–594) survive. Use Harmony / negative-population normalization. |
| Per-donor median univariate gating | INS fraction-positive varies 8× across donors; Beta recall drops 80% donor-to-donor. |
| Per-cell QC of flagged cells | Hundreds of thousands of cells → thousands of hours. QC must be **cluster-level**. |

## 3. The core empirical finding — this panel's gates are unreliable in BOTH directions

Per-(marker, donor) 2-component GMM gating (the basis of `napariGater`'s seeds and of scimap's
auto-gating) is **not trustworthy on this panel**, measured directly on donor 6374:

- **Over-gated** (called positive in ~90–100 % of *all* cells): `SOX2` (90 %), `EpCAM`, `Ker8_18`,
  `Vimentin`. Because `Ductal = anypos(Keratin_5, TP63, SOX2)`, the bad `SOX2` gate alone dragged
  **72 % of cells into "Ductal"** — biologically impossible (pancreas is acinar-dominated).
- **Under-gated** (called positive in only ~1–5 % even of clearly-epithelial clusters):
  `Pan_Cytokeratin`. Since `Acinar` requires `PanCK pos`, acinar cells are *lost*.
- **Well-gated** (calibrated to the true ~1 % positive mode): `INS, GCG, SST, Keratin_5, TP63, CD3e,
  CD8, CD4, FOXP3, CD68, CD163, CD31, CD34, …` — the low-prevalence, well-separated markers.

This is exactly why manual napari gating was painful: no univariate threshold (Otsu / GMM valley /
crossover) can know `SOX2` should be rare or `PanCK` common without a biological prior or multivariate
context. **No purely univariate method will be robust here.**

## 4. Method comparison on donor 6374 (200k cells)

| Method | Script | Ductal call | Result | Verdict |
|---|---|---|---|---|
| Per-donor GMM gating | `auto_gate_per_donor.py` | **72 %** (wrong) | Mechanically correct, faithful rescale (99.999 % vs scimap), but SOX2 over-gate poisons it | Fast unblock only; needs ~15 curated gate overrides |
| Leiden + hierarchical annotation | `cluster_annotate.py` | **16 %** (plausible); Acinar 24 % (correct) | **Groups cells correctly** (real endocrine/endothelial/stromal clusters); annotation is the hard, semi-supervised part | **Recommended foundation** |
| Astir (probabilistic, no gates) | `run_astir.py` | — | Reference-free, GPU-fast; at 30 **and** 150 epochs collapsed to ~46 % Unknown + 31 % "Other" + 18 % Smooth muscle (SMA/Vimentin high-background dominate after standardization) | Defeated out-of-the-box by the **same marker-quality issue** as gating; needs marker curation (drop high-background channels) before it is usable |

**Key result:** multivariate clustering *fixed the SOX2 problem the moment it grouped cells* — a
SOX2-high-but-otherwise-acinar cell clusters with acinar, so that cluster's mean SOX2 stays low and it
is correctly called Acinar. The clusters themselves are real biology; **what needs human input is
labelling ~40–80 clusters on a heatmap, not 59 sliders.**

## 5. Recommended workflow (staged)

**Phase 0 — fast unblock (no napari):** `auto_gate_per_donor.py` computes per-(marker, donor) gates,
rescales per donor, runs the existing scimap rules, and emits a **per-(marker, donor) diagnostic PDF**.
QC = scan the PDF, correct the ~10–15 mis-gated markers (e.g. `SOX2`, `PanCK`) in a small overrides CSV
(`--gate-overrides`), re-run. Produces `obs['phenotype']` immediately for downstream `aggregate_islets.py`.

**Phase 1 — robust multivariate typing (the accuracy core):**
1. **Cluster** (`cluster_annotate.py`): arcsinh → per-marker z → PCA → **Harmony(imageid)** (batch
   correction; `harmonypy` already installed) → Leiden. Reliable, batch-robust grouping.
2. **Auto-annotate** clusters by **dominant-lineage** scoring of cluster-mean z against the rules
   hierarchy, with `Unassigned (low/ambiguous)` for non-specific/artifact clusters and a top-1-vs-top-2
   `annot_margin` confidence. Outputs a cluster→phenotype table + a cluster×marker **heatmap** (the QC
   tool) + per-cell labels.
3. **Astir cross-check** (`run_astir.py`, `rules_to_astir_yaml.py`): reference-free per-cell posteriors
   with no gating; reuses the rules as a marker YAML. Independent second opinion + per-cell confidence.
4. **Consensus + cluster-level QC:** where clustering and Astir (and curated-gate Phase 0) agree → high
   confidence; disagreements and `Unassigned`/low-margin clusters → the ~40–80 cluster-level decisions a
   human reviews on the heatmap.

**Phase 2 — validation without ground truth:** marker coherence (Beta=INS⁺, T=CD3e⁺…), spatial
homotypic coherence (`X/Y_centroid`), per-donor composition trends (β-loss + immune-gain ND→Aab+→T1D),
cross-method concordance (~60–80 % realistic), and an expert spot-check of a few islets per donor.

**Phase 3 — downstream:** run `aggregate_islets.py` unchanged (`obs['phenotype']` contract preserved).

## 6. Scaling to the full 23 M cells

Cluster a stratified per-donor subsample (≈150k/donor ≈ 2.25 M) with Harmony+Leiden, annotate the
clusters once, then propagate to all 23 M cells with a fast classifier on the normalized markers
(no RAPIDS needed). Astir minibatches over the full set on GPU. Phase-0 gating + rescale are vectorised
(the fast rescale is validated faithful) and run on the full set directly.

## 7. How to run

```bash
# Phase 0 — fast per-donor gating (scimap env)
conda run -n scimap python scripts/senior/auto_gate_per_donor.py \
    --out data/singlecell_protein_phenotyped.h5ad \
    --gates-out data/per_donor_gates.csv \
    --diagnostics-pdf data/per_donor_gate_diagnostics.pdf
# (review the PDF; correct mis-gated markers in overrides.csv; re-run with --gate-overrides overrides.csv)

# Phase 1a — cluster + annotate (scvi-env)
conda run -n scvi-env python scripts/senior/cluster_annotate.py \
    --max-cells-per-donor 150000 --harmony --resolution 2.0 \
    --out-prefix data/cluster_subsample

# Phase 1b — Astir cross-check (astir-env)
conda run -n scvi-env python scripts/senior/rules_to_astir_yaml.py --out scripts/senior/astir_markers.yaml
conda run -n astir-env python scripts/senior/run_astir.py --epochs 200 \
    --markers scripts/senior/astir_markers.yaml --out-prefix data/astir_full
```

Envs: `scimap` (scimap 2.3.5, scanpy, harmonypy) · `scvi-env` (scanpy 1.11, leidenalg, harmonypy,
torch+CUDA) · `astir-env` (astir 0.1.5, torch+CUDA, pyarrow).

> **⚠️ CORRECTION (2026-06-17, user ground truth):** Sections 8b–8d below describe a workflow that was *built and run mechanically*, but these data have **NOT** undergone successful, quality-controlled phenotyping. There was **no genuine expert cluster annotation** — `cluster_labels_curated.csv` was not a real expert review, and the "13/13 marker-coherent / validated end-to-end" results are a weak no-ground-truth coherence check, **not** sign-off. Treat every "completed / curated / validated" claim in §8b–8d as **draft/automated only**; the current app phenotypes are provisional. The methodology and empirical findings in §1–§7 remain valid. Genuine QC'd phenotyping (clustering + a real expert review; immune subtyping via Junior's scimap gating method; Astir skipped) is the approved plan `~/.claude/plans/review-the-recent-plans-sprightly-lagoon.md`.

## 8b. Marker-reliability audit + multi-donor clustering (draft — NOT expert-QC'd)

`marker_reliability_audit.py` scored all 59 markers across the 15 donors with two signals: gate
statistics (per-donor GMM positive-fraction, separation, cross-donor CV) and a **gate-free
cluster-informativeness** signal (`top_cluster_z` = how strongly the marker concentrates in some
cluster). The latter is decisive — it separates "bad antibody" from "legitimately common":

- **Confirmed background → drop from features + signatures** (`AUDIT_BACKGROUND` in `cluster_annotate.py`,
  `--exclude-markers` for clustering, `--exclude-markers SOX2 …` for the Astir YAML):
  **SOX2** (53 % positive, top_cluster_z 1.02), **EpCAM**, **E_cadherin**, **DAPI**, plus housekeeping
  **Beta_actin, Bcl_2, PCNA, Ki67**.
- **Excellent markers that are merely *under-gated*** (low positive-fraction, very high cluster
  specificity — proof clustering beats gating): INS 3.45, CD3e 3.28, GCG 2.90, HLA_DR 2.90, MPO 2.81,
  B3TUBB 2.86, CD31 2.57, Keratin_5 2.59, SST 2.48, CD68 2.02.
- **Weak immune subtype markers** (CD4 1.28, FOXP3 1.23, CD20 1.43) — immune subtyping is the hardest.

Output: `marker_audit_marker_audit.csv` + `_distributions.pdf` (per-marker per-donor densities).

**Curated 15-donor clustering** (`cluster_annotate.py --harmony --exclude-markers <8 background>`,
750k cells, Harmony bug fixed): with the final annotation (dominant-lineage best-guess that always
descends to a leaf by *own* markers, with a lineage-strength floor + a near-tie → Unassigned rule, and
audit markers dropped from signatures) the cluster-level call distribution is biologically sensible —
**Acinar 27 %, Ductal 13 %, Blood vessel 6 %, M2 macrophage 6 %, Fibroblast 4 %, Beta 3 %, Neural 3 %,
CD8 T 1 %, …**, with **~36 % flagged `Unassigned`** (low-signal / near-tie clusters) for heatmap review.
Each cluster also gets **subtype alternatives** (e.g. `Beta cell(3.41); Alpha cell(2.65); Delta cell(1.27)`)
so the human review is a ranked choice, not a blank. This is the intended cluster-level QC: ~40 cluster
decisions on the heatmap, not 59 sliders × 15 donors.

**Final hand-off loop:** review `cluster_subsample_cluster_annotation.csv` against
`cluster_subsample_heatmap.png` → save a curated `leiden,phenotype` CSV → `finalize_phenotype.py`
trains a classifier on the labelled subsample and propagates to all 23 M cells (+ fuses Phase-0 gate /
Astir labels into a confidence + flag) → `validate_phenotypes.py` for the no-ground-truth checks →
`aggregate_islets.py` unchanged.

## 8d. Draft run with auto-curated cluster labels (NOT a genuine expert review)

The 40 clusters were read against the full panel (high-specificity audit markers) and relabelled in
`cluster_labels_curated.csv` (+ `_with_rationale.csv`) — fixing several auto-calls (cluster 30 = clearly
endothelial; 13 & 25 = Granulocyte by CD66/MPO; 3 = myeloid by CD206/CD209/CD68; 8 = Ductal by TP63),
with 14 genuinely low-signal / proliferating / debris clusters kept `Unassigned` for image review.
`finalize_phenotype.py --cluster-labels cluster_labels_curated.csv` propagated to all **23.28 M cells**
(mean confidence **0.872**) -> `singlecell_protein_phenotyped.h5ad`.

**Endocrine subtype split.** Resolution-2.0 clusters lump islet endocrine cells together (clusters 28/31
are INS⁺GCG⁺SST⁺), so a single cluster label cannot separate Beta/Alpha/Delta. Since INS/GCG/SST are the
panel's most reliable markers (audit top_cluster_z 3.4/2.9/2.5), `finalize_phenotype.py` splits the
endocrine cells **per-cell by dominant hormone** (INS=Beta, GCG=Alpha, SST=Delta) — 864,386 endocrine
cells → **Beta 432,714 / Alpha 370,984 / Delta 60,688** (realistic islet composition; alpha elevated as
expected in a T1D-weighted cohort).

Validation (`validate_phenotypes.py`): **13/13 phenotypes marker-coherent** — Beta=INS +3.15, Alpha=GCG
+2.54, Delta=SST +2.11, plus all non-endocrine types. The split **unmasks the central T1D biology** that
was blurred when β/α/δ were merged: **β-cell fraction ND 3.1 % → Aab+ 1.4 % → T1D 0.3 %** (≈10× monotonic
loss — near-complete beta destruction in T1D) and **immune fraction ND 9.3 % → Aab+ 4.8 % → T1D 18.4 %**
(insulitis). Remaining refinements: resolve the rare immune subtypes (CD4 T, B, NK, DC) folded into larger
clusters (higher resolution or targeted immune subclustering), and review the 14 `Unassigned` clusters on
the image.

## 8c. Full 23 M scale-up — ran mechanically (NOT validated by expert QC)

`finalize_phenotype.py` trained a classifier on the 484k labelled subsample cells and propagated to all
**23,280,387 cells** in chunks (mean classifier confidence **0.906**; 13 GB phenotyped h5ad written),
proving the clustering path scales on one box. The *draft* (pre-curation) over-extends Ductal/Acinar
because the 36 % Unassigned subsample clusters get force-assigned — re-running after the heatmap cluster
review (curated `leiden,phenotype` CSV) is what produces the final phenotyping.

## 8. Bottom line

The 59-plex panel cannot be phenotyped by any plug-and-play automated method because several markers are
intrinsically poorly gated (SOX2/EpCAM over; PanCK/INS/CD3e under) — triangulated across gating, Astir,
and the marker audit. The robust path (scripts implemented; the human cluster-review step has **not** genuinely been done — see §8b correction) is: **cluster reliably
(multivariate, Harmony-batch-corrected, audit-curated features) → auto-suggest cluster labels + subtype
alternatives + confidence → review ~40 clusters on the heatmap → propagate to 23 M + fuse methods +
validate.** This replaces the 59-slider × 15-donor napari grind with ~40 heatmap-guided cluster
decisions. All seven scripts are built, compile, and have been run on real data; the one irreducibly
human step is the cluster-label review (domain expertise), which the pipeline is built to consume and
re-run.
