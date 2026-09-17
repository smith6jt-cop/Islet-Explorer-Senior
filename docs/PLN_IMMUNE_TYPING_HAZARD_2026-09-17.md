# pLN immune typing — the gDNA hazard, measured per section

**The ruling is settled.** 2026-09-17, the maintainer: *"All immune cell calls from pancreas can be
applied to pLN."* The immune lineage vocabulary established on pancreas — `T_NK` / `Immune_T`,
`B_Plasma` / `Immune_B`, `Myeloid`, `Granulocyte` — is the same vocabulary in the pLN, and the typing
method that produced the pancreas calls is licensed to run on pLN cells. This document does not
reopen that. It exists so that a pLN immune number, once it exists, is read with the right error
bar beside it.

Everything below was measured on 2026-09-17 from files on this machine. Nothing is quoted from a
summary.

---

## 1. What the ruling licenses, and the two boundaries it does not move

Licensed:

- The same lineage labels, with the same meanings, in both tissues.
- Running the typer on pLN cells.

Not licensed, and not extended here:

- **No spatial claim crosses the tissue boundary.** The pancreas and pLN sections come from
  different blocks mounted on one slide. They share a pixel grid and share nothing anatomically.
  No distance, adjacency, neighbourhood, KDE or shared scatter may span them. A lymph node has no
  islet, so a signed distance to an islet there is **unavailable**, not zero.
- **No label rows are copied.** Measured on 6442: its 958,744 pancreas calls match 100% of its
  `panc` cells and 44 of its 177,799 `pln_1` cells, and those same 44 are exactly the `cell_id`
  collisions between the two bundles. The two bundles were segmented separately, so the `cell_id`
  namespaces are disjoint and the 44 are coincidence, not correspondence. The typer must be **run**
  on the pLN cells.
- **Immune only.** Endocrine, Exocrine and the other non-immune lineages are outside the ruling.

---

## 2. The model, reproduced before it was applied

`docs/XENIUM_PIPELINE.md` §0 and §2 give the vendor model and the resulting core-size bias. It was
re-derived here from the constants in that document, and it reproduces every published row before
being carried to any new section.

**The model.**

1. Every gene probe picks up genomic DNA at the same rate a `genomic_control` probe does. That rate
   is `λ` = `genomic_control_probe_counts_per_control_per_cell`, read from the bundle's
   `metrics_summary.csv`. `λ` is the **expected gDNA counts per gene per cell**.
2. Under the detection predicate now in force — `raw count >= 1`, hardcoded at six places in
   `analysis/xenium/scripts/hierarchical_phenotype.py` — a single gene is spuriously detected by
   gDNA alone with probability `p = 1 − exp(−λ)` (Poisson).
3. A lineage's identity core is anchored when `ANCHOR_MIN_GENES = 2` of its `n` core genes are
   detected. gDNA is uniform across probes, so the number spuriously detected is
   `Binomial(n, p)`, and

   ```
   P(core spuriously anchored by gDNA alone) = P( Binomial(n, p) >= 2 )
   ```

**Reproduction against `XENIUM_PIPELINE.md`.** Using that document's own λ for the flagged pair
(0080870, λ = 0.0252021):

| quantity | model here | document |
|---|---:|---:|
| P(≥1 count from gDNA, one gene) | 0.02489 | 2.49% |
| P(≥2 counts from gDNA, one gene) | 0.000312 | 0.031% |
| P(≥1) for 0081366 (λ = 0.02571) | 0.02538 | 2.54% |
| P(anchor) Exocrine, n=8 | 0.0157 | 0.016 |
| P(anchor) Endocrine, n=18 | 0.0728 | 0.073 |
| P(anchor) Immune_B, n=58 | 0.4250 | 0.425 |
| P(anchor) Vascular, n=87 | 0.6405 | 0.641 |
| P(anchor) Myeloid, n=111 | 0.7663 | 0.766 |
| P(anchor) Immune_T, n=129 | 0.8337 | 0.834 |
| Immune_T under a 2-counts-per-gene floor | 0.000784 | 0.0008 |

Every row reproduces. The model is the document's, not an invention.

**Core sizes, recovered from `core/panels`** (the one taxonomy authority), and they reproduce the
document's sizes exactly:

| lineage as the ruling names it | files unioned | n core genes |
|---|---|---:|
| `Immune_T` / `T_NK` | `03_immune_T_cell.csv` ∪ `05_immune_NK_ILC.csv` (101 + 31, 3 shared) | **129** |
| `Myeloid` | `06_immune_myeloid.csv` ∪ `07_immune_granulocyte.csv` (90 + 21, 0 shared) | **111** |
| `Immune_B` / `B_Plasma` | `04_immune_B_plasma.csv` | **58** |
| `Granulocyte` alone | `07_immune_granulocyte.csv` | **21** |
| Endocrine (reference) | `tissues/pancreas/identity_core/01_pancreas_endocrine.csv` | 18 |
| Exocrine (reference) | `tissues/pancreas/identity_core/02_pancreas_exocrine.csv` | 8 |

The Xenium release's `label_space.json` folds `Granulocyte` into `Myeloid`; the ruling names it
separately, so it is carried in the table as its own 21-gene core as well.

**Panel, read from this cohort's own bundles, not assumed.** `gene_panel.json` in the pLN bundles:
**5,101 genes, 21 `genomic_control`, 40 `negative_control`, 3,291 `deprecated_codeword`**. This is
5,101 genes, not the 5,001 of `XENIUM_PIPELINE.md` §6 — the cohort carries a 100-gene `hMulti_100g`
custom add-on on top of `hAtlas_v1.1`. `n_probes = 21` is unchanged. Checked on three bundles
(0041139, 0041140, 0059947); the h5ads carry 5,101 `var` rows, consistent.

**Vendor reproduction, this cohort.** `λ × 5,101` against the vendor's
`estimated_number_of_false_positive_transcripts_per_cell_including_genomic_counts`, all 46 bundles:
**maximum relative error 0.0000%**. Validation check 1 of `XENIUM_PIPELINE.md` §5 passes on every
bundle in this cohort, pLN included.

---

## 3. The h5ads cannot supply this number — the bundle metrics must

Trap 1 of `XENIUM_PIPELINE.md` holds for every pLN input. Checked on all 24
`/data/xenium/h5ad_out/*_pLNXenium.h5ad`:

```
n_var = 5101 in all 24;  Intergenic_Region_* = 0;  NegControl* = 0;  Deprecated* = 0
```

`spatialdata_io.xenium` loaded Gene Expression only, so the genomic controls are **not in the
h5ads at all**. There is no way to measure a gDNA rate from these files, and any attempt to infer
one from gene counts would be measuring expression. The rate comes from `metrics_summary.csv` or it
is unavailable.

**Coverage.** `/data/xenium_pln/_metrics/` holds 46 `metrics_summary.csv`. Joining them to
`core/phenocycler/data/integration/xenium_paths.csv` (serial + region → donor + roi, the repo's
authority) resolves **22 pLN sections**. Two of the 46 do not resolve to any declared bundle:
`0059865__PanT` and `0059907__PanT`. Those are the two serials whose declared pLN bundles —
`0059865__pLN__20250328__202338` (6539) and `0059907__pLN__20250328__202338` (6436) — are **not on
this machine in any form**: not in `/data/xenium_pln/`, not in `/data/xenium_bundles/`, not in
`/data/xenium/processed/`. Only `Panc` and `PanT` were pulled for those two serials.

So: **22 of the 24 pLN h5ads have a measured rate; 6436 and 6539 are UNMEASURED.** They are
unmeasured, not zero and not low — nothing here says anything about their background.

A cell-count correspondence links each h5ad to its bundle (h5ad `n_obs` sits just under the
bundle's `num_cells_detected`, QC having removed the difference) for all 22. For the two
unmeasured ones the nearest bundle by cell count is the `PanT` section of the same serial
(6436: 370,852 against 373,800; 6539: 218,842 against 222,443), which is suggestive of a naming
question upstream but is **not** a measurement of their gDNA and is not used as one below.

---

## 4. Per-section table

`λ` is the quantity the model needs. `×0041323` is λ against the sample every detection constant in
this repo was tuned on (λ = 0.000325). `adj. rate` is
`adjusted_genomic_control_probe_rate`, carried for continuity with
`docs/XENIUM_GDNA_TRIAGE_2026-09-17.md`. `λ pLN / λ panc` is against that donor's own pancreas
section on the same serial. `P(...)` is `P(Binomial(n_core, 1 − exp(−λ)) ≥ 2)` — the per-cell
probability that the lineage's identity core is anchored **by genomic DNA alone**, under the
`raw count >= 1` predicate now in the code.

| pLN h5ad | donor | roi | serial | λ (gDNA counts/gene/cell) | ×0041323 | adj. rate | λ pLN / λ panc | median tx | bundle local | P(Immune_T) | P(Myeloid) | P(B_Plasma) | P(Granulocyte) |
|---|---|---|---|---:|---:|---:|---:|---:|---|---:|---:|---:|---:|
| `6347` | 6347 | pln_1 | 0041140 | 0.002492 | 7.67× | 0.05094 | 5.40× | 217 | yes | **0.0415** | 0.0316 | 0.00934 | 0.001261 |
| `6580` | 6580 | pln_1 | 0041139 | 0.001037 | 3.19× | 0.05522 | 6.55× | 69 | yes | **0.0081** | 0.0061 | 0.00171 | 0.000223 |
| `6457a` | 6457 | pln_1 | 0037637 | 0.000901 | 2.77× | 0.05475 | 4.86× | 49 | metrics only | **0.0062** | 0.0046 | 0.00130 | 0.000168 |
| `6374` | 6374 | pln_1 | 0041032 | 0.000781 | 2.40× | 0.01707 | 0.63× | 200 | yes | **0.0047** | 0.0035 | 0.00098 | 0.000127 |
| `6579a` | 6579 | pln_1 | 0041035 | 0.000775 | 2.39× | 0.01865 | 2.19× | 193 | metrics only | **0.0046** | 0.0035 | 0.00096 | 0.000125 |
| `6482` | 6482 | pln_1 | 0037354 | 0.000556 | 1.71× | 0.01217 | 2.26× | 203 | yes | **0.0024** | 0.0018 | 0.00050 | 0.000064 |
| `6442` | 6442 | pln_1 | 0059947 | 0.000454 | 1.40× | 0.00467 | 0.86× | 559 | yes | **0.0016** | 0.0012 | 0.00034 | 0.000043 |
| `6476` | 6476 | pln_1 | 0037634 | 0.000417 | 1.28× | 0.00891 | 2.95× | 217 | metrics only | **0.0014** | 0.0010 | 0.00028 | 0.000036 |
| `6579b` | 6579 | pln_2 | 0041035 | 0.000378 | 1.16× | 0.00957 | 1.06× | 176 | metrics only | **0.0011** | 0.0008 | 0.00023 | 0.000030 |
| `6461` | 6461 | pln_1 | 0059847 | 0.000351 | 1.08× | 0.01100 | 0.69× | 150 | yes | **0.0010** | 0.0007 | 0.00020 | 0.000026 |
| `6493` | 6493 | pln_1 | 0037191 | 0.000335 | 1.03× | 0.00799 | 2.06× | 190 | yes | **0.0009** | 0.0007 | 0.00018 | 0.000024 |
| `6457b` | 6457 | pln_2 | 0037637 | 0.000317 | 0.98× | 0.00926 | 1.71× | 154 | metrics only | **0.0008** | 0.0006 | 0.00016 | 0.000021 |
| `6579c` | 6579 | pln_3 | 0041035 | 0.000305 | 0.94× | 0.01011 | 0.86× | 132 | metrics only | **0.0007** | 0.0006 | 0.00015 | 0.000019 |
| `6380` | 6380 | pln_1 | 0037502 | 0.000293 | 0.90× | 0.01118 | 1.32× | 99 | metrics only | **0.0007** | 0.0005 | 0.00014 | 0.000018 |
| `6551` | 6551 | pln_1 | 0037376 | 0.000291 | 0.90× | 0.00355 | 3.19× | 432 | yes | **0.0007** | 0.0005 | 0.00014 | 0.000018 |
| `6533` | 6533 | pln_1 | 0041393 | 0.000283 | 0.87× | 0.01078 | 2.46× | 105 | yes | **0.0006** | 0.0005 | 0.00013 | 0.000017 |
| `6523` | 6523 | pln_1 | 0037499 | 0.000283 | 0.87× | 0.01048 | 1.55× | 133 | metrics only | **0.0006** | 0.0005 | 0.00013 | 0.000017 |
| `6521` | 6521 | pln_1 | 0037199 | 0.000278 | 0.85× | 0.00674 | 1.09× | 184 | yes | **0.0006** | 0.0005 | 0.00013 | 0.000016 |
| `6593` | 6593 | pln_1 | 0059945 | 0.000224 | 0.69× | 0.00289 | 1.30× | 387 | yes | **0.0004** | 0.0003 | 0.00008 | 0.000010 |
| `HDL115` | 115 | pln_1 | 0041117 | 0.000213 | 0.66× | 0.00753 | 1.24× | 150 | yes | **0.0004** | 0.0003 | 0.00007 | 0.000009 |
| `6566` | 6566 | pln_1 | 0059968 | 0.000197 | 0.61× | 0.01065 | 1.93× | 78 | yes | **0.0003** | 0.0002 | 0.00006 | 0.000008 |
| `HDL112` | 112 | pln_1 | 0041396 | 0.000190 | 0.58× | 0.00504 | 1.05× | 168 | yes | **0.0003** | 0.0002 | 0.00006 | 0.000008 |
| `6436` | 6436 | pln_1 | 0059907 | **UNMEASURED** | — | — | — | — | absent | **UNMEASURED** | — | — | — |
| `6539` | 6539 | pln_1 | 0059865 | **UNMEASURED** | — | — | — | — | absent | **UNMEASURED** | — | — | — |

Reference points, same model, for scale:

| section | λ | P(Immune_T) | P(Myeloid) | P(B_Plasma) | P(Exocrine) |
|---|---:|---:|---:|---:|---:|
| 0041323 — what every constant was tuned on | 0.000325 | 0.0008 | 0.0006 | 0.0002 | 0.0000 |
| **worst measured pLN (6347)** | **0.002492** | **0.0415** | **0.0316** | **0.0093** | **0.0002** |
| 0080870 panc (6533, flagged Pheno-Xeno pair) | 0.025202 | 0.8337 | 0.7663 | 0.4250 | 0.0157 |
| 0081366 panc (112, flagged Pheno-Xeno pair) | 0.025710 | 0.8419 | 0.7758 | 0.4349 | 0.0163 |

---

## 5. Three things the measurement changes

**a. The 0.834 headline does not transfer to any pLN section in this cohort.** It belongs to
λ = 0.0252. The worst measured pLN λ is 0.002492 — **9.9× lower** — and the worst measured pLN
P(Immune_T) is 0.0415, twenty times below the flagged pair. The composition argued in the task
brief ("largest cores × worst gDNA") is directionally right and is visible in the table: the
immune lineages remain the most exposed at every λ, and the pLN is the noisier tissue. But the
magnitude is an order of magnitude smaller than the flagged pair, and it was worth measuring
rather than carrying the headline across.

**b. `adjusted_genomic_control_probe_rate` is not the quantity this model needs, and it reorders
the sections.** The adjusted rate is λ divided by the section's own per-gene depth, so a shallow
section is pushed up the ranking by its shallowness. Under λ:

- 6580 pLN is **first** on the adjusted rate (0.05522, 12.2× tuned) and **second** on λ (3.19×).
- 6347 pLN is **third** on the adjusted rate (0.05094, 11.2×) and **first** on λ (7.67×).
- 6551 pLN sits at 0.00355 adjusted (0.90× on λ) yet carries 3.19× its own pancreas λ.

Spurious anchoring depends on absolute background per gene per cell, which is λ. Shallow depth is
a separate problem — it shrinks real signal against a fixed background — and it is why a
group-level call is the right architecture, but it does not itself raise P(anchor). Any partition
drawn off the adjusted rate is drawn off the wrong axis.

**c. The pLN-is-noisier finding holds on λ too, at a smaller multiple.** On λ, the pLN section is
noisier than its own pancreas section on the same serial in **18 of 22** paired sections, median
**1.63×**, max **6.55×** (6580). `XENIUM_GDNA_TRIAGE_2026-09-17.md` reported 17 of 18, median 3.9×,
max 39× on the adjusted rate. Both are true of their own quantity; the λ version is the one the
anchoring model consumes.

Two bookkeeping notes on that triage, for whoever owns it: it counts 23 pLN rows against the 22
that resolve here, because it classifies `0041140__Panc_pLN` as `pln` while `xenium_paths.csv`
declares that bundle as 6347's **panc**; and `cassette_name` in `metrics_summary.csv` is not a
donor id — 0041393 carries `nPOD6553` where the path authority says 6533, and 0059907 and 0059945
both carry `6436`. Neither affects the pLN rows in the table above, which are keyed on
serial + region through `xenium_paths.csv`.

---

## 6. The partition

The pass condition is not invented here. `XENIUM_PIPELINE.md` §5 check 3 states it:
**P(Immune_T satisfied by gDNA) → < 0.01**, Exocrine approximately unchanged. Applied to the
measured λ of each pLN section:

**Group A — existing floors are fine (21 of 22 measured).**
`6580`, `6457a`, `6374`, `6579a`, `6482`, `6442`, `6476`, `6579b`, `6461`, `6493`, `6457b`,
`6579c`, `6380`, `6551`, `6533`, `6523`, `6521`, `6593`, `HDL115`, `6566`, `HDL112`.
All carry P(Immune_T) < 0.01 under `raw count >= 1`; the largest is 6580 at 0.0081, which is at
81% of the bar and is the one to watch. Fifteen of these twenty-one sit below 0.0015.

**Group B — run `gdna.py` first (1 of 22 measured).**
`6347` (`pln_1`, serial 0041140). P(Immune_T) = **0.0415**, 4.2× the bar; P(Myeloid) = 0.0316;
P(B_Plasma) = 0.0093. At λ = 0.002492 its expected gDNA background is 12.71 counts per cell
against a median of 217 transcripts — 5.9% of a median cell. Its bundle
(`/data/xenium_pln/output-XETG00298__0041140__pLN__20241025__204655`) is on this machine, so
Stage 1 can run on it as specified in `XENIUM_PIPELINE.md` §2.

**Group C — UNMEASURED, and blocked on a file, not on a decision (2 of 24).**
`6436` and `6539`. Their declared pLN bundles are not on this machine. Until a
`metrics_summary.csv` for `0059907__pLN` and `0059865__pLN` is pulled, their λ is unavailable and
they cannot be assigned to A or B. Do not default them into A.

A caution about Group A: "existing floors are fine" means only that **this one hazard** clears
**this one bar**. It is not a clearance to type. The other blockers in `XENIUM_PIPELINE.md` are
untouched by anything here — see §7.

Headroom. Solving the model for the bar, **P(Immune_T) = 0.01 at λ = 0.0011561**. Because
P(≥2) ≈ C(n,2)·λ² at small λ, a factor-of-two error in λ moves P(Immune_T) by roughly a factor of
four, so the margins are worth stating as λ ratios rather than as P differences:

| section | λ | λ(bar) / λ |
|---|---:|---:|
| 6580 | 0.001037 | 1.11× |
| 6457a | 0.000901 | 1.28× |
| 6374 | 0.000781 | 1.48× |
| 6579a | 0.000775 | 1.49× |
| 6482 | 0.000556 | 2.08× |
| 6442 | 0.000454 | 2.54× |
| 6476 | 0.000417 | 2.77× |
| the remaining 14 | ≤ 0.000378 | ≥ 3.06× |

6580 crosses the bar if its λ is understated by 11%.

---

## 7. What a pLN immune number will and will not support

**It will support**, once the typer has been run on pLN cells:

- A per-donor, per-section count of cells called `T_NK`, `B_Plasma`, `Myeloid` or `Granulocyte`
  within one pLN section, reported as a raw number with its denominator.
- A within-donor comparison of that donor's pancreas and pLN immune composition, stated as two
  separate measurements from two separate blocks — never as one distribution and never with a
  distance between them.
- A cross-donor comparison **among pLN sections of the same group**, with the λ column carried
  beside it, and with 6347 either corrected or reported separately.

**It will not support:**

- Any spatial statement spanning pancreas and pLN. Distance to an islet in a pLN is unavailable.
- Pooling a donor's multiple pLN sections into one number without saying so. 6457 has two
  (`pln_1`, `pln_2`) and 6579 has three (`pln_1` … `pln_3`); they are separate sections with
  separate λ, differing here by 2.8× within 6457 and 2.5× within 6579.
- A cross-donor immune-fraction comparison that mixes Group A with uncorrected 6347: at
  P(Immune_T) = 0.0415 against a Group A median near 0.0008, 6347's immune fraction carries a
  different spurious floor from everyone else's, and a difference that follows λ is a measurement
  of λ.
- **Any check of the form "the pLN came out mostly immune, so the typing worked."** A lymph node
  being mostly lymphocytes is the expectation, not evidence. A result that passes by reproducing
  it has measured the expectation.

**Cohort bookkeeping that applies to the table.** 6457 is centrally excluded
(`core/phenocycler/phenocycler/cohort.py` — technical, poor DAPI), so `6457a` and `6457b` are
measured here but are not live: **20 live measured pLN sections**. 6579 was reinstated on
2026-09-17 in that same file; its pancreatitis is not retracted, only its exclusion, so `6579a/b/c`
must be stratified or reported separately rather than pooled. 6591 declares no `pln_1` at all —
recorded unavailable, not missing.

**Blockers that this document does not clear and that no λ value can clear:**

1. **Validation check 2 is unrun for every section here.** `XENIUM_PIPELINE.md` requires any
   detection predicate to put zero-template Y probes at or below ~0.3%. That check needs the full
   bundle, needs the donor's sex, and needs the control probes — which are absent from every
   h5ad. Eight of the 22 measured pLN sections have no local bundle at all
   (`6523`, `6380`, `6476`, `6457a`, `6457b`, `6579a`, `6579b`, `6579c`), so the check cannot be
   run on them from this machine.
2. **The scale problem is untouched.** `XENIUM_PIPELINE.md` measures `raw count >= 1` at 1.88%
   false detection on a zero-template probe against 0.55% for `lognorm >= log1p(1)`, and states the
   fix is a change of scale rather than a higher floor. The repo's raw-count floors are the
   vulnerable ones and are still in place.
3. **A 2-counts-per-gene floor collapses this hazard entirely, everywhere.** At 6347's λ — the
   worst measured — P(Immune_T anchored by gDNA) falls from 0.0415 to **7.9 × 10⁻⁸**. The hazard
   quantified in this document is specific to the `>= 1` predicate.
4. **The `nullgate` null already absorbs part of this.** `match_strata` matches on per-gene
   detection rate and mean, and gDNA is uniform per probe, so a matched random core carries the
   same background as the real core and it largely cancels in the contrast. The numbers in §4 are
   the hazard to a **per-cell** anchoring decision, which is what
   `hierarchical_phenotype.py` makes today. A group-level call gated against a detection-matched
   null is exposed to much less of it — and recording the background-predicted contrast per group,
   as §2 Stage 4 of the pipeline document specifies, would make how much explicit.

One observation, dated: `/data/xenium/processed/0041393/pln_1/0041393_phenotyped.h5ad` was written
at 14:07 on 2026-09-17, during this session. Until today `/data/xenium/processed/<serial>/` held
only `panc/`. That file is 6533's pLN — Group A, P(Immune_T) = 0.0006. Nothing here inspects it; it
is noted so that its provenance is not mistaken for something older.

---

## 8. Reproducing

Environment: `/home/smith6jt/IO60panc2nd/Islet-Explorer-Senior/.venv/bin/python` (pandas, scipy,
h5py). No GPU, no `scvi-env`, no `xenium_analysis` needed — nothing here loads a model.

Inputs, all read-only:

| what | where |
|---|---|
| λ per bundle | `/data/xenium_pln/_metrics/*.metrics_summary.csv`, field `genomic_control_probe_counts_per_control_per_cell` |
| donor + roi per bundle | `core/phenocycler/data/integration/xenium_paths.csv`, keyed on the final path component |
| identity core sizes | `core/panels/data/identity_core/`, `core/panels/data/tissues/pancreas/identity_core/` |
| panel composition | `gene_panel.json` in any `/data/xenium_pln/output-XETG00298__*` bundle, `type.descriptor` counts |
| control absence in the h5ads | `/data/xenium/h5ad_out/*_pLNXenium.h5ad`, `var/_index` |

The computation, in full:

```python
import math, pandas as pd
from scipy.stats import binom

lam = float(pd.read_csv(metrics_csv).iloc[0]["genomic_control_probe_counts_per_control_per_cell"])
p   = 1 - math.exp(-lam)                 # P(one gene detected by gDNA alone, raw count >= 1)
P   = 1 - binom.cdf(1, n_core, p)        # P(core anchored by gDNA alone, ANCHOR_MIN_GENES = 2)
```

with `n_core` = 129 (`Immune_T`/`T_NK`), 111 (`Myeloid` incl. granulocyte), 58
(`Immune_B`/`B_Plasma`), 21 (`Granulocyte` alone). Constants: `n_genes = 5101`, `n_probes = 21`,
λ(0041323) = 0.000325, λ(0080870) = 0.0252021, λ(0081366) = 0.02571.

Nothing was written outside this document.
