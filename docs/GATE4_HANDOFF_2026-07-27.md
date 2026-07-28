# Handoff — 2026-07-27 (evening): Gate 4 complete on 20 donors

> **SUPERSEDED — go to [`GATE5_HANDOFF_2026-07-28.md`](GATE5_HANDOFF_2026-07-28.md).** That handoff
> carries the current state, the one blocking decision, and the corrections to this document.


> **SUPERSEDED IN PART by [`GATE5_ROUND_2026-07-27.md`](GATE5_ROUND_2026-07-27.md) — read that next.**
> The Gate-5 evidence package is built (`data/gate5/`), the per-image review bundle is built
> (`expanded_review_v12/`, 34 PDFs + 161 overlays), and Gate 6's cell accounting is proven exact.
> Four corrections to this document:
> 1. **The islet-core Immune disease gradient is aliased with acquisition BATCH** (Kruskal p=0.0016 by
>    batch vs 0.0094 by status; batch 2 holds every high-immune donor and zero ND donors). The Fusion
>    version was a perfect proxy but is **not** the mechanism. Open item 7's "immune environment as
>    disease progresses" cannot be claimed from this cohort without a batch-aware analysis.
> 2. **CD11b's positives are 81.1% epithelial-marker⁺ and 80.8% epithelial⁺ outside islets**, with
>    78.6% carrying no other immune marker. An earlier claim in this banner that CD11b is
>    "under-called ~4×" is **RETRACTED** — there is no ground truth for the pilot's 7–10%, and raising
>    CD11b would add more of the same cells. Its divisor is still the only one that ROSE under
>    `reference_only`/p99 (1.46–2.25×), which remains unexplained.
> 3. **Open item 4's SST framing is withdrawn** — it is a per-donor pass/fail (bimodal `target_fold`,
>    nothing between 1.81 and 2.54), not a disease effect (Fisher p=0.167).
> 4. **Open item 5 should cover 6380, 6534, 6623 as well as 6566** — INS's arm separation collapses on
>    all four, 6380 worst of all, yet they get three different dispositions.

**Start here.** Supersedes the Gate-4 sections of
[`ANATOMY_CHECK_2026-07-27.md`](ANATOMY_CHECK_2026-07-27.md) (still authoritative for the anatomy
evidence and the divisor sweep) and `HANDOFF_2026-07-27_qc_and_v10.md` (still authoritative for cell
QC, REDSEA v10 and the duplicate-detection finding — note it is **gitignored**, so it exists only on
this machine).

Commits: parent through `1def2d2`, submodule `redsea-restore-curated-5compartment`. 242 tests.

---

## Where things stand

**Gate 4 is complete and clean on all 20 eligible donors.** 480 donor-pairs: **479 VALID**, 1
`NO_TARGET_EXPRESSION_CONFIRMED` (6566 INS), **zero technical failures, 0.0% Unresolved**.
25,151,028 cells typed.

| region | Epithelial | Endothelial | Immune | Mesenchymal | Neural | `Other` |
|---|---|---|---|---|---|---|
| islet core | **78.9%** | 11.2% | 8.8% | 0.4% | 0.1% | **0.7%** |
| peri-islet | 58.3% | 10.6% | 19.7% | 1.9% | 0.4% | 9.0% |
| tissue | 48.3% | 13.3% | 21.7% | 2.0% | 0.5% | 14.2% |

Islet-core Epithelial is **60.5–85.5% across all 20 donors** (median 83.1%), against **1–6%** under
the 7-marker gates. Cohort `Other` is **13.7%**, from ~65% in the rejected pilot.

Islet-core cells are endocrine by QuPath annotation and hormone-validated (core INS 407 / GCG 1004
against 10 / 71 in tissue, cohort medians). Those annotations are drawn from the image and are
independent of every marker threshold, so this is the one external check the pipeline has.

## What changed to get here

1. **Gate expansion, 7 → 24 markers** (maintainer sign-off, `pair_reviews_accepted.csv`). `Other` was
   the ordered residual's terminal bucket, so a cell whose only lineage marker was ungated could not
   be typed — B cells, DC, neutrophils, macrophage subsets, endocrine cells, lymphatics and smooth
   muscle were landing there by construction.
   - Immune: CD3e CD68 CD11b **CD79a CD11c Iba1 CD206 MPO CD4 CD8**
   - Epithelial: E_cadherin **Pan_Cytokeratin Ker8_18 EpCAM Keratin_5 INS GCG SST**
   - Endothelial: CD31 **CD34 Podoplanin** · Neural: B3TUBB · Mesenchymal: Vimentin **SMA**
2. **Policy frozen**: `reference_only` / `p99`. Interior choice — `max` calls E-cadherin on 4.3–18.9%
   and cannot run at all; `p95` puts 7–31% of islet-core cells in Immune.
3. **`COMPARTMENT_ORDER` → Immune > Endothelial > Epithelial > Neural > Mesenchymal.** Adding the
   hormones to Epithelial inverted a validated control: 67.8% of islet-core CD31⁺ cells were being
   claimed by Epithelial, and 78% of them were also hormone⁺ (against 1.8% in tissue). Islet
   capillaries sit in a mass of hormone-bright neighbours. Both original constraints still hold.
4. **Feature compartments resolve per PAIR, not per run.** 14 targets have only whole-cell values, so
   a run-level `Cell` setting had been degrading the 10 pairs that do have Membrane. That alone was
   the CD68/6523 failure.

## Open items, in priority order

1. **The gold standard is the only thing that can verify any of this.** Bundle is built and blinded at
   `data/gold_standard_v1/` — 2,500 cells, 5 donors spanning all disease groups and acquisition
   batches, 3 weighted strata, QuPath import/export groovy, scorer that reports inter-annotator κ
   first. **Needs annotator time, not code.** Protocol: `scripts/gold_standard/README.md`.
2. **All 14 new markers reached production on screen diagnostics alone.** The per-image biological
   review is outstanding for every one. `SHORTLIST_PAIRS` (34) carries a rationale per pair naming
   what to look for. Build the bundle from `screen_level1_expansion`.
3. **Tissue Endothelial 13.3%**, up from 2–6% under CD31 alone. **CD34 is the likely cause** — called
   at 16.8% cohort-median and expressed on fibroblasts as well as endothelium. Flagged, not
   reordered around.
4. **SST anomaly, unexplained.** target_fold median 1.68 in ND with low SBR on 62% of ND donors,
   against 20.7 in T1D and 55.2 in AAB. Delta cells should not be weakest in non-diabetic pancreas.
5. **INS on 6566 is `CONFIRMED_ABSENT` provisionally**, from screen diagnostics not image review. INS
   target_fold collapses 199/137 (AAB/ND) → 8.4 (T1D) with low SBR on 44% of T1D donors — beta-cell
   loss. Confirm from the image before reporting any INS-dependent result.
6. **Diagnostic flags across the run**: 18 `target_low_sbr`, 104 `reference_undersized` of 480 rows.
   Reported, never gating.
7. **Gates 5–7**: 20-donor dashboard + per-donor sign-off, then broad-lineage acceptance, then rebuild
   descendants (H5AD, islet aggregation, trajectory, app).

## Things a future session will otherwise get wrong

- **The proportion rule is inadmissible for majority-positive targets.** Pan-CK, Ker8_18, EpCAM *and*
  E-cadherin are undersized against every reference on 100% of donors (ratios 0.08–0.18). That is a
  property of the targets, not the references — a marker positive on 25–50% of cells cannot have a
  larger negative control. The `E_cadherin ← CD31` exception is the general case, not a special one.
  `reference_fold` is 19–37× for all of them; only the *relative* size criterion fails.
- **The QuPath channel normalization does not apply to RESTORE.** `results/norm_stats.csv`,
  `norm = (raw − p_low)/span`, six channels. RESTORE fits per image, so pure scaling is provably inert
  (÷294, ÷100, ÷1000 all reproduce the same result exactly). Only the offset does anything, by
  perturbing the binned Otsu floor, and the answer swings 42% on its value. It also puts 37% of cells
  negative into a pipeline whose NNMF needs nonnegative input. The columns never reached any parquet.
- **`restore_apply --on-technical-failure`**: `halt` (default) stops the run; `unavailable` records
  the donor-marker as blocking so its cells become `Unresolved`. Neither absorbs a failure as a
  biological negative. Currently no failures, so the default is fine again.
- **`figures/*` and `HANDOFF*` are gitignored.** The anatomy-check PNGs and the QC handoff are
  local-only. Every number they support is in the tracked markdown.
- **`data/cells/` is heterogeneous**: 22 partitions, but only the 20 eligible donors are QC-filtered.
  An ad-hoc `read_parquet('data/cells/donor_id=*')` mixes pre- and post-QC cells.

## Reproduce

```bash
conda activate scvi-env
python -m pytest Phenocycler_Analysis/tests -q                       # 242
python -m phenocycler.restore_apply --config config.ini --jobs 8      # ~2 h, 24 pairs x 20 donors
python -m phenocycler.lineage --config config.ini --jobs 8
# evidence
python scripts/anatomy_check/run_arms.py && python scripts/anatomy_check/figures.py
Phenocycler_Analysis/scripts/run_divisor_sweep.sh <out> --pair-set <set> --arm reference_only/p99
python -m phenocycler.reference_selection --results <screen> --output <dossier>
```

Key artifacts: `data/restore_gated_redsea/` (tri-state calls + `apply_manifest.csv`),
`data/phenotype/broad/` (compartments), `data/restore_pair_validation/reference_selection_v2/`,
`screen_v11_qc/`, `screen_level1_expansion/`, `sweep_qc/`, `data/gold_standard_v1/`,
`data/donor_acquisition_batch.csv`.
