# pLN Xenium — immune typing by cluster-then-name, as run

**Maintainer's method, 2026-09-17, in four statements.** *"All immune cell calls from pancreas can
be applied to pLN."* *"Separate out the immune cells with canonical markers/genes and then use
subpanels."* *"Clustering works well in pLN because of the regional cell types."* *"Tissue does not
mean anything at all."*

The last one is the load-bearing one, and it is literal. In the pancreas a QuPath pixel classifier
hands each cell a context — islet, blood vessel, lymphatic, nerve, ordinary tissue — and a reference
level is fitted against it. A lymph node section is one block of node: there is no partition of the
section, so there is no "ordinary tissue" to draw a reference population from and the context
machinery has nothing to attach to. **Clustering is the substitute for context.** The partition is
discovered rather than supplied, and each cluster is then named by which canonical genes stand out
*against the other clusters on that same section*.

| | |
|---|---|
| Producer | `analysis/xenium/scripts/pln_immune_typing_v2.py` (new) |
| Release | `data/development/xenium_typing/release_pln_immune_v2/` (new; the draft is left in place) |
| `release_id` | `xenium-pln-immune-cluster-named-dev-20260917`, `schema_version` 1, `roi` `pln_1` |
| Clusters consumed from | `data/development/typing/pln/clusters/xenium/assignments/` |
| Cells | **6,007,789**, 22 sections, **20 donors** |
| Latch | development; nothing here is a release of record |

`data/development/xenium_typing/release/` — the pancreas release — was not touched, and
`release_pln_immune/` — the draft this supersedes — was not touched either.

---

## 1. What is withdrawn, and is not in this build

1. **Any number derived from pancreas cells reaching a pLN call** — level, threshold, scale factor,
   quantile. The producer reads no pancreas data of any kind. The earlier attempt that transferred
   each donor's pancreas p99 onto their pLN slide is withdrawn; its numbers are a property of that
   mistake and are not cited here as findings.
2. **`subpanel.fit_levels` / `reference_level` / `context == 'tissue'` on the pLN.** There is no
   context column and no reference population, and none is synthesised.
3. **`IMMUNE_REQUIRE_MIN = 0.10` as an absolute floor on pLN values.** It is defined against a
   reference level that does not exist here.
4. **The draft release's detection predicate**, `raw >= depth / section_median_depth`. That is a
   raw-count floor that rises with cell depth, applied in the tissue where gDNA background is worst.
   §4 says what replaced it and shows the replacement binds nowhere.

---

## 2. Inputs, measured

**The partition is consumed, never recomputed.** `core/phenocycler/scripts/pln/cluster_pln.py`
partitioned each section within donor; this lane reads its `assignments/` verbatim and re-derives
only the per-cluster gene profile, because it needs its own detection predicate (the clustering
track's `det__` columns are `lognorm > 0` under a section-median normalisation, which is not the
predicate declared here). Cluster ids are the clustering track's own `<section>_cNN`.

- **22 sections, 20 donors, 364 clusters.** 6579 contributes three sections (`6579a/b/c`); every
  other donor one. 6457 was refused upstream by `phenocycler.cohort.filter_eligible_donors` and no
  lane here restates that list.
- **6,007,789 cells in the h5ads**; 5,920,568 carry a cluster. The remainder is 45,171 cells below
  the clustering track's QC (`nCount_Xenium < 10` or `nFeature_Xenium < 5`) and 42,050 excluded for
  having no nucleus (§4).
- **429 genes wanted, 429 present, 0 absent, on every one of the 22 sections** — the five immune
  identity cores, the four non-immune cores, the stage-1 canonicals, the Y probes and the regional
  set, all read from `core/panels/data/identity_core/`. No copy of those CSVs was written.
- **`CD3D` and `LYZ` are absent from all 24 pLN h5ads.** Nothing is substituted for them; the
  release records them as unavailable.
- **No control feature of any kind exists in the h5ad `var` index** — no `NegControl*`,
  `Intergenic*`, `Unassigned*` or `Deprecated*` — measured on all 22 sections. `spatialdata_io`
  loads Gene Expression only, so **no control probe was ever there to exclude**, and nothing in this
  release claims controls were excluded.
- `X` row sums equal `obs['nCount_Xenium']` for **all 6,007,789 cells**, so depth is the cell's own
  measured total and not a second quantity wearing the same name.

**Not in this release, and why.** `6450`, `6505`, `6534`, `6591`, `6623` have no pLN h5ad on this
machine. `6457` is centrally excluded. Nothing about those donors is reported as zero.

---

## 3. The method, as run

### Stage 0 — what "stands out" means

For gene *g* and cluster *k* of one section, with `mean(k,g)` the cluster's mean lognorm and
`det(k,g)` its detection fraction:

```
standout(k,g)  <=>  z(k,g) >= 1      where z = (mean(k,g) - mean over that section's clusters) / SD over them
                AND det(k,g) > det(rest of the section, g)
```

Both halves are contrasts **inside one section**. Nothing is imported: no level, no quantile, no
threshold that would have to mean the same thing on two slides. The z half asks whether the cluster
is high for that gene relative to the section's own clusters; the detection half stops a difference
in a handful of high-count cells from carrying a gene on its own.

### Stage 1 — immune against non-immune, on clusters

The immune set is the maintainer's list: `PTPRC CD3E CD3G CD4 CD8A MS4A1 CD79A CD19 CD68 ITGAX
ITGAM NKG7 KLRD1`. The structural set is the naming genes of the four non-immune identity cores
(`11_endothelial_vascular`, `12_pericyte_smooth_muscle`, `15_neural_glial`, `16_epithelial_general`
— 52 on panel).

```
Immune             <=>  >= 2 immune canonicals stand out  AND  f_imm > f_str
Non_immune         <=>  >= 2 structural naming genes stand out  AND  f_str > f_imm
otherwise          ->   Uncalled_ambiguous, with the reason recorded
```

`f` is the **fraction** of each set standing out, not the count: the sets are 13 and 52 genes and a
count comparison would hand the larger one an advantage that is not biology.

| stage 1 | clusters | cells | median cluster depth | median immune standouts | median structural standouts |
|---|---:|---:|---:|---:|---:|
| Immune | 137 | 3,072,781 | 201 | 4 | 1 |
| Non_immune | 165 | 2,125,949 | 194 | 0 | 12 |
| Uncalled_ambiguous | 62 | 721,838 | **71** | 0 | 1 |

The 62 ambiguous clusters are the low-depth ones — cell-weighted median depth 120 against 229 for
the called clusters — and 44 of the 62 reach two standouts on neither side. They are named, not
smoothed away.

### Stage 2 — identity within immune, from the subpanels

Sources: `core/panels/data/identity_core/03_immune_T_cell.csv` (101 genes, 20 of which may name),
`04_immune_B_plasma.csv` (58 / 26), `05_immune_NK_ILC.csv` (31 / 22), `06_immune_myeloid.csv`
(90 / 31), `07_immune_granulocyte.csv` (21 / 19). Which genes **may name** a lineage is read off the
panel's own `evidence` column — `canonical` or `both` may name, `cellxgene` ("marks lineage (CL)")
scores and describes but cannot admit a label.

Score of a core = mean `z` over that core's on-panel genes; the winner is the argmax.

> **The small-core defect is the thing to fix**, and the fix is not a margin. A core with few genes
> can win an argmax on near-zero evidence, so the winner is admitted only when **at least one gene
> that may name that lineage actually stands out in that cluster**. A cluster that fails is
> `Immune_uncalled`, named, and **never reassigned to the runner-up** — absent evidence is not
> evidence for something else.

The gate fires on **1 cluster of 137, 17,390 cells** (`6476_c11`, winning core `Immune_B_Plasma`,
zero standout naming genes). It fires rarely because a cluster mean is a far less noisy statistic
than a single cell's counts — which is the whole reason the draft's per-cell argmax needed a floor
of 2 and still produced the calls verification refuted. **136 of 136 called clusters carry ≥ 1
standout naming gene of their own lineage** (minimum by lineage: Myeloid 3, Granulocyte 5, NK_ILC 7,
T 7, B_Plasma 11).

The evidence is auditable gene by gene. Every Granulocyte cluster, with the naming genes that stood
out in it:

| donor | cluster | cells | standout naming genes of the granulocyte core |
|---|---|---:|---|
| 6461 | `6461_c07` | 4,124 | CCR3 CEACAM3 CSF3R CXCR1 CXCR2 FCGR3B MMP8 MPO PRTN3 RNASE3 |
| 6521 | `6521_c13` | 4,671 | CEACAM3 CTSG CXCR1 CXCR2 ENPP3 GATA2 IL5RA MPO PRG2 |
| 6551 | `6551_c09` | 2,540 | CCR3 CEACAM3 CSF3R CXCR1 CXCR2 FCGR3B MMP8 MPO RNASE3 |
| 6579 | `6579c_c06` | 4,091 | CCR3 CEACAM3 CSF3R CXCR1 CXCR2 MMP8 MMP9 MPO PRTN3 RNASE3 |
| 6580 | `6580_c11` | 5,927 | CEACAM3 CSF3R CXCR1 CXCR2 FCGR3B |
| 6593 | `6593_c18` | 6,158 | CCR3 CEACAM3 CSF3R CXCR1 CXCR2 FCGR3B MMP8 PRG2 |

`6521_c13` stands out on `GATA2`, `IL5RA`, `PRG2`, `ENPP3` and `CTSG` and on none of the neutrophil
genes the other five carry. That is the eosinophil/mast branch of the granulocyte core, not
neutrophils, and the label does not split it. The 17 NK_ILC clusters carry 7–17 naming genes each
(`KLRD1`, `KLRF1`, `KLRC1/2`, `PRF1`, `NCR1`, `EOMES`, `TBX21`, `KIR*`); two of them —
`6593_c15` and part of `6374_c01` — lean on `GATA3`/`KIT`/`RORC`/`IL1RL1`, which is the ILC branch
of that core rather than NK. `Immune_NK_ILC` is the union and is not split here.

### Inheritance and heterogeneity

Cells inherit their cluster's call verbatim. What that costs is measured per cluster and written
into `clusters/cluster_calls.parquet`:

| call | clusters | cells | cells whose own core argmax agrees | cells detecting ≥1 naming gene of the call |
|---|---:|---:|---:|---:|
| Immune_T | 45 | 1,421,211 | 74.2% | 91.4% |
| Immune_B_Plasma | 33 | 983,313 | 93.8% | 97.6% |
| Myeloid | 35 | 414,866 | 64.6% | 88.8% |
| Immune_NK_ILC | 17 | 208,490 | **12.0%** | 62.9% |
| Granulocyte | 6 | 27,511 | **34.7%** | 57.9% |

Per-cell argmax is the draft's own estimator, and it is reported here beside the cluster call rather
than in place of it. Where the two disagree most — NK_ILC and Granulocyte — the disagreement is
visible in the release, not hidden behind a majority label.

---

## 4. Detection, the nucleus rule, and the controls

### Detection

```
lognorm(c,g) = log1p(raw(c,g) * 1e4 / depth(c)),   detected  <=>  lognorm >= log1p(1)
```

The target is **fixed at 1e4**, not the section median. For any cell of depth ≤ 10,000 the predicate
is exactly `raw >= 1`, so no depth-scaled floor operates. Measured: **0 of 6,007,789 cells have
depth above 10,000**, and **0 of 119,234,523 non-zero gene entries are dropped** by the target. The
predicate is `raw >= 1` on this data, and that is measured rather than asserted.

### A cell with no nucleus is excluded from every calculation

The h5ads carry no nucleus geometry, so the only handle is the vendor bundle `cells.parquet`, matched
by donor **and** by subset containment of the h5ad's own `cell_id` set — never by filename order.

- **14 of 22 sections have a local bundle**: 42,050 cells with `nucleus_count == 0` are excluded
  from the cluster aggregates as well as from the calls, and are labelled `Unavailable_no_nucleus`.
- **8 sections have none** — `6436`, `6476`, `6493`, `6523`, `6539`, `6579a/b/c`. For those the rule
  is **UNMEASURABLE**, recorded as such in `release.json`, never as passed. Transferring those
  bundles' `cells.parquet` is the whole of what it takes to close that gap.

### The Y-probe control

`TBL1Y`, `USP9Y` and `EIF1AY` are the floor group: zero-template in a female donor and, unlike
`DDX3Y` and `KDM5D`, not cross-hybridising with a highly expressed X paralog. In a female donor
every detection of them is false, so their rate is the background this predicate admits. **Nine
female donors, eleven sections** (the brief said two; measured, it is nine).

| section | donor | floor group | all five Y probes |
|---|---|---:|---:|
| HDL115 | 115 | **0.0029%** | 0.028% |
| 6579b | 6579 | 0.0026% | 0.030% |
| 6533 | 6533 | 0.0052% | 0.029% |
| 6476 | 6476 | 0.0062% | 0.050% |
| 6523 | 6523 | 0.0064% | 0.030% |
| 6579c | 6579 | 0.0065% | 0.040% |
| HDL112 | 112 | 0.0069% | 0.038% |
| 6380 | 6380 | 0.0073% | 0.027% |
| 6374 | 6374 | 0.0148% | 0.096% |
| 6579a | 6579 | 0.0172% | 0.088% |
| 6482 | 6482 | **0.0173%** | 0.084% |

**0.0026–0.0173% against the ~0.3% requirement — inside it by a factor of 17 to 114.** In the eleven
male sections the same probes run 3.49–23.17%, which is template, not background, and is why they
are not a control there: for male donors this quantity is **UNMEASURED**, and the release says so
rather than writing 0.

### gDNA

`data/development/xenium_gdna/gdna_triage_2026-09-17.parquet` carries a measured pLN rate for 18 of
the 20 donors — 0.0029 (6593) to 0.0552 (6580); `6347` (0.0509) and `6580` (0.0552) are the two
above 10× the tuned rate. **`6436` and `6539` are UNMEASURED** and are recorded as the string
`UNMEASURED`, never as 0.

---

## 5. The release

```
data/development/xenium_typing/release_pln_immune_v2/
├── release.json                       schema_version 1, release_id, roi pln_1, donor list,
│                                      the method block, and sha256 for all 44 inputs
│                                      (22 h5ads + 22 cluster assignment parquets)
├── label_space.json                   10 labels, each with the reason it exists
├── calls/donor_id=<d>/roi=pln_1/data.parquet
├── clusters/cluster_calls.parquet|csv 364 clusters: calls, standouts, scores, heterogeneity
├── clusters/cluster_gene_stats.parquet 156,156 rows: per cluster per gene mean, z, detection,
│                                      detection in the rest of the section, standout flag
├── per_donor_summary.csv · controls.csv · section_reports.json
├── verification_narrow_canonicals.csv the adversarial re-check of §7
└── compare_draft_summary.csv · compare_draft_donor_<d>.csv · compare_draft.json
```

`release.json`'s method block describes what the code does, statement by statement — the draft's was
refuted for claiming things its code did not do, so every clause there is one a reader can check
against `pln_immune_typing_v2.py`.

**Keying.** Rows are keyed on `cell_id` within `donor_id` / `roi`, and `cell_id` is the Xenium
barcode verbatim, `-1` suffix included. **One measured exception: 6579's three sections share 35
`cell_id` values**, so for that donor the key is (`section_key`, `cell_id`). `release.json` records
this at `duplicate_cell_ids_within_donor`; no other donor has any.

### Per donor, every label including the uncalled

| donor | status | sex | cells | clusters | Immune | % | Non_immune | Ambiguous | not clustered | no nucleus | Immune_T | Immune_B_Plasma | Immune_NK_ILC | Myeloid | Granulocyte | Immune_uncalled |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 112 | ND | F | 221,578 | 18 | 116,262 | 52.5 | 65,242 | 37,175 | 506 | 2,393 | 54,428 | 44,768 | 2,853 | 14,213 | 0 | 0 |
| 115 | ND | F | 247,399 | 13 | 146,065 | 59.0 | 53,289 | 44,294 | 1,287 | 2,464 | 65,381 | 49,540 | 15,074 | 16,070 | 0 | 0 |
| 6347 | AAb+ | M | 79,522 | 16 | 33,118 | 41.6 | 33,970 | 11,664 | 160 | 610 | 11,435 | 14,155 | 1,161 | 6,367 | 0 | 0 |
| 6374 | ND | F | 164,364 | 16 | 67,142 | 40.8 | 76,671 | 18,844 | 496 | 1,211 | 32,743 | 19,148 | 7,669 | 7,582 | 0 | 0 |
| 6380 | T1D | F | 264,010 | 11 | 142,654 | 54.0 | 73,975 | 42,879 | 2,050 | 2,452 | 73,479 | 47,325 | 0 | 21,850 | 0 | 0 |
| 6436 | T1D | M | 370,852 | 12 | 223,226 | 60.2 | 147,364 | 0 | 262 | n/a | 100,728 | 42,963 | 42,296 | 37,239 | 0 | 0 |
| 6442 | ND | M | 176,137 | 22 | 106,988 | 60.7 | 40,133 | 25,105 | 105 | 3,806 | 50,764 | 35,089 | 3,252 | 17,883 | 0 | 0 |
| 6461 | ND | M | 109,476 | 17 | 51,504 | 47.0 | 53,758 | 2,882 | 615 | 717 | 13,941 | 21,265 | 5,428 | 6,746 | 4,124 | 0 |
| 6476 | ND | F | 290,429 | 16 | 154,752 | 53.3 | 113,405 | 21,719 | 553 | n/a | 47,823 | 72,382 | 6,527 | 10,630 | 0 | **17,390** |
| 6482 | ND | F | 168,551 | 15 | 81,303 | 48.2 | 62,040 | 24,045 | 474 | 689 | 39,256 | 15,285 | 10,223 | 16,539 | 0 | 0 |
| 6493 | ND | M | 161,178 | 14 | 75,015 | 46.5 | 30,075 | 55,590 | 498 | n/a | 42,019 | 26,909 | 0 | 6,087 | 0 | 0 |
| 6521 | AAb+ | M | 479,006 | 17 | 194,343 | 40.6 | 232,140 | 47,639 | 1,513 | 3,371 | 75,545 | 61,065 | 15,478 | 37,584 | 4,671 | 0 |
| 6523 | T1D | F | 600,995 | 12 | 340,260 | 56.6 | 209,966 | 44,109 | 6,660 | n/a | 172,929 | 125,858 | 0 | 41,473 | 0 | 0 |
| 6533 | T1D | F | 97,237 | 13 | 43,865 | 45.1 | 37,703 | 13,932 | 892 | 845 | 22,395 | 16,871 | 0 | 4,599 | 0 | 0 |
| 6539 | ND | M | 218,842 | 17 | 116,062 | 53.0 | 71,292 | 31,047 | 441 | n/a | 46,506 | 37,103 | 10,937 | 21,516 | 0 | 0 |
| 6551 | T1D | M | 174,918 | 15 | 69,062 | 39.5 | 104,620 | 0 | 315 | 921 | 42,017 | 15,820 | 0 | 8,685 | 2,540 | 0 |
| 6566 | T1D | M | 692,878 | 26 | 335,686 | 48.4 | 191,650 | 135,854 | 20,218 | 9,470 | 198,200 | 97,693 | 18,336 | 21,457 | 0 | 0 |
| 6579 | T1D | F | 577,283 | 36 | 267,112 | 46.3 | 227,125 | 78,504 | 4,542 | n/a | 121,880 | 79,740 | 11,722 | 49,679 | 4,091 | 0 |
| 6580 | T1D | M | 119,600 | 39 | 38,454 | 32.2 | 51,494 | 25,393 | 3,098 | 1,161 | 12,747 | 14,748 | 245 | 4,787 | 5,927 | 0 |
| 6593 | T1D | M | 793,534 | 19 | 469,908 | 59.2 | 250,037 | 61,163 | 486 | 11,940 | 196,995 | 145,586 | 57,289 | 63,880 | 6,158 | 0 |
| **20** | | | **6,007,789** | **364** | **3,072,781** | | **2,125,949** | **721,838** | **45,171** | **42,050** | **1,421,211** | **983,313** | **208,490** | **414,866** | **27,511** | **17,390** |

`n/a` in the nucleus column is the eight sections with no local bundle — unmeasurable, not zero.
Stage 1 and stage 2 each sum to 6,007,789 exactly.

---

## 6. Against the draft release

The draft (`release_pln_immune/`, 13 donors, 3,569,071 bundle cells) is left in place and is
superseded rather than overwritten. On the cells both cover — the 13 shared donors, h5ad cells only
— rows are the draft's `lineage_native`, columns this release's `lineage`:

| draft ↓ / v2 → | Immune_T | Immune_B_Plasma | Immune_NK_ILC | Myeloid | Granulocyte | Non_immune | Ambiguous | not clustered | no nucleus |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Immune_T (614,484) | 413,912 | 3,014 | 59,916 | 19,497 | 108 | 68,088 | 49,928 | 21 | 0 |
| Immune_B_Plasma (593,442) | 25,884 | 439,387 | 4,417 | 10,162 | 387 | 55,188 | 57,998 | 19 | 0 |
| Immune_NK_ILC (14,040) | 4,479 | 197 | 5,886 | 471 | 245 | 1,873 | 889 | 0 | 0 |
| Myeloid (106,195) | 2,141 | 968 | 2,248 | 81,517 | 2,520 | 14,134 | 2,659 | 8 | 0 |
| Granulocyte (4,708) | 287 | 71 | 221 | 928 | 1,874 | 1,085 | 241 | 1 | 0 |
| Immune_unresolved (623,846) | 250,268 | 38,390 | 44,353 | 59,503 | 5,496 | 118,511 | 106,922 | 403 | 0 |
| Non_immune (788,959) | 46,343 | 25,722 | 8,035 | 29,348 | 3,405 | 644,501 | 31,583 | 22 | 0 |
| Ambiguous (738,928) | 72,533 | 43,284 | 11,932 | 24,966 | 9,385 | 349,367 | 197,770 | 29,691 | 0 |
| Unavailable (39,598) | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 39,598 |

What this says, and only this:

- **The two lineages verification indicted move the most.** Of the draft's 4,708 Granulocyte calls,
  1,874 (39.8%) remain Granulocyte; of its 14,040 NK_ILC calls, 5,886 (41.9%) remain. Neither is
  reassigned wholesale to one place — they spread across T, Myeloid, Non_immune and Ambiguous.
- **The draft's unresolved and ambiguous cells are where most of the new calls come from.** Its
  623,846 `Immune_unresolved` and 738,928 `Ambiguous` supply 250,268 + 72,533 T calls and
  38,390 + 43,284 B calls here. A cluster mean survives sparsity that a single cell's counts do not;
  that is the mechanism, and it cuts both ways (§7).
- **The no-nucleus set agrees exactly** — the draft's 39,598 `Unavailable` cells inside the h5ad map
  one to one onto `Unavailable_no_nucleus`.
- Stage 1 agrees on **2,123,248 of the 3,524,200 shared cells (60.2%)**; stage 2 on **1,587,077
  (45.0%)**. Agreement is not itself evidence for either release — it is stated so the size of the
  change is on the record.

Coverage differs too: this release types **20 donors** because it works from the h5ad and the
cluster assignment, where the draft needed the vendor bundle for `nucleus_count` and so stopped at
13. The trade is stated in §4 — seven of those seven extra donors have no measurable nucleus rule.

---

## 7. What these numbers do not say

**The per-cell evidence problem is relocated, not solved, and here is the measurement.** The check
that refuted the draft, recomputed from the h5ads against this release: for each cell, does it detect
any of the narrow canonical list for the lineage it was called?

| call | cells | detecting **none** of its lineage's narrow canonicals |
|---|---:|---:|
| Immune_T | 1,421,211 | 61,984 (**4.4%**) |
| Immune_B_Plasma | 983,313 | 53,399 (**5.4%**) |
| Myeloid | 414,866 | 73,916 (**17.8%**) |
| Immune_NK_ILC | 208,490 | 107,040 (**51.3%**) |
| Granulocyte | 27,511 | 16,326 (**59.3%**) |

For NK_ILC and Granulocyte that is the same rate the draft was refuted at (50.5% and 55.4%). The
difference is what the claim attaches to: the draft called each cell individually, so a cell with no
canonical of its own lineage was a call resting on nothing. Here the call is a property of the
**cluster** — every one of those clusters stands out on 5 to 17 naming genes of its lineage against
the section's other clusters — and the per-cell rate is the cost of inheritance, reported in the
table above, in `verification_narrow_canonicals.csv`, and per cluster in `cluster_calls.parquet`. A
reader who needs cell-level evidence should filter on `n_naming_detected_for_cluster_call`, which is
written per cell for exactly that purpose. **The granulocyte and NK_ILC counts here are cluster-level
statements. They should not be read as 27,511 individually corroborated granulocytes.**

**Thin stage-1 calls are flagged, not hidden.** 23 Immune clusters (239,643 cells) and 22
Non_immune clusters (293,754 cells) rest on exactly two standouts — the minimum the rule admits.
`stage1_thin_two_standouts` marks them and `standout_immune_genes` / `standout_structural_genes`
name the genes, so a reader can see that e.g. a Non_immune cluster carried by `CDH1` and `PRPH`
alone is not the same evidence as one carried by seventeen endothelial genes.

**The standout threshold is a choice, and its effect is measured.** Re-deciding everything at
z ≥ 0.75 and z ≥ 1.5 (recorded per section in `section_reports.json`):

| label | z ≥ 0.75 | **z ≥ 1.0 (released)** | z ≥ 1.5 |
|---|---:|---:|---:|
| Immune_T | 1,465,872 | 1,421,211 | 1,054,357 |
| Immune_B_Plasma | 986,948 | 983,313 | 978,039 |
| Immune_NK_ILC | 211,685 | 208,490 | 192,555 |
| Myeloid | 414,866 | 414,866 | 385,018 |
| Granulocyte | 37,215 | 27,511 | 29,582 |
| Non_immune | 2,193,423 | 2,125,949 | 1,852,656 |
| Uncalled_ambiguous | 591,196 | 721,838 | 1,428,361 |
| Immune_uncalled | 19,363 | 17,390 | 0 |

B_Plasma and Myeloid barely move; T and the ambiguous pool trade against each other; Granulocyte is
the small-count label and is the least stable of the five.

**Internal coherence, reported; expectation, not.** A lymph node is expected to be mostly immune and
**that expectation is not a check** — nothing here was accepted because it reproduced it, and 6580
at 32.2% immune and 6521 at 40.6% are reported unchanged beside 6442 at 60.7%. What *is* checkable is
contradiction. Cell-weighted detection rates by called lineage:

| called | PTPRC | CD3E | CD4 | CD8A | MS4A1 | CD79A | CD19 | CD68 | ITGAX | NKG7 | MPO |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Immune_T | 70.1 | 57.8 | 19.1 | 22.1 | 8.8 | 9.9 | 4.5 | 3.7 | 0.7 | 3.1 | 0.05 |
| Immune_B_Plasma | 66.6 | 17.3 | 7.3 | 5.7 | **77.7** | **61.0** | **59.5** | 4.7 | 1.8 | 1.5 | 0.06 |
| Myeloid | 61.1 | 24.1 | 44.3 | 8.9 | 9.6 | 9.2 | 5.7 | **36.1** | **12.9** | 2.9 | 0.19 |
| Immune_NK_ILC | 79.1 | 59.5 | 14.0 | **54.0** | 11.7 | 10.5 | 5.5 | 5.7 | 2.5 | **23.6** | 0.06 |
| Granulocyte | 52.1 | 6.7 | 7.6 | 2.5 | 4.8 | 4.7 | 3.3 | 6.5 | 19.2 | 1.3 | 0.20 |

**CD4 19.1% and CD8A 22.1% within the T-called cells — 41.2% together, no contradiction.** The
withdrawn attempt put CD4 at 69% and CD8 at 58% on the same cells, which was one. Two things the
table shows that are not contradictions but are worth naming: `CD4` is higher in Myeloid-called
cells (44.3%) than in T-called ones, which is what CD4 does on monocytes and macrophages; and the
NK_ILC-called cells carry CD3E at 59.5% and CD8A at 54%, which says these clusters sit against the
cytotoxic T compartment and that the union label is doing work there.

---

## 8. What the method does not license

- **No spatial claim crosses the pancreas/pLN boundary.** They are different blocks on one slide.
  Nothing here computes a distance, adjacency, neighbourhood, KDE or shared scatter, and
  `obsm['spatial']` is absent from all 24 h5ads in any case. The pancreas/pLN comparison this feeds
  is **compositional**, joined on `donor_id` and disease status and nothing else.
- **No non-immune lineage is resolved.** `Non_immune` is terminal. The four non-immune cores appear
  only as far as stage 1 needs them to exclude.
- **Nothing claims acinar or ductal.** Exocrine is one label in this project, and it would mean
  nothing in a node.
- **The regional genes are evidence, not calls.** `CR2`, `CXCL13`, `CCL19`, `CCL21`, `BCL6`,
  `AICDA`, `SELL`, `CD27`, `MKI67` are carried in `cluster_gene_stats.parquet` per cluster. No
  follicle, germinal-centre or T-zone call is made from them here; the regional reading is
  `docs/PLN_REGIONAL_CLUSTERING_2026-09-17.md`'s.
- **No cluster gated anything.** The partition describes; it admits no label on its own and refuses
  none.

---

## 9. Reproduce

```bash
cd /home/smith6jt/IO60panc2nd/Islet-Explorer-Senior
.venv/bin/python analysis/xenium/scripts/pln_immune_typing_v2.py --stage call
.venv/bin/python analysis/xenium/scripts/pln_immune_typing_v2.py --stage compare
# one section or donor:  --only 6533
```

CPU, repo `.venv`, no scVI and no GPU. **~2 minutes wall** for all 22 sections including re-reading
every h5ad and hashing all 44 inputs. The producer takes no random draw of any kind; re-running
`--only 6533` reproduced its `data.parquet` **byte for byte** (sha256 identical). Note that a
`--only` run rewrites `release.json`, the cluster tables and the summaries for that subset alone —
re-run the full `--stage call` afterwards, which is what was done here.

## 10. Files

| path | |
|---|---|
| `analysis/xenium/scripts/pln_immune_typing_v2.py` | **the producer** |
| `data/development/xenium_typing/release_pln_immune_v2/**` | the release |
| `docs/PLN_IMMUNE_XENIUM_2026-09-17.md` | this file |
| `data/development/xenium_typing/release_pln_immune/**` | the draft — **left in place**, superseded |
| `data/development/xenium_typing/release/**` | the pancreas release — untouched |
| `core/panels/data/identity_core/*.csv` | read only; no copy written |
