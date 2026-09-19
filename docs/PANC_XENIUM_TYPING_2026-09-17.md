# Pancreas Xenium 5K typing — the scoring defect, the fix, and the cohort release

Everything below was measured on 2026-09-17 on this machine by running something. Nothing is
quoted from a summary and no number is carried in from another tissue, another donor or another
lane.

Owner files: `analysis/xenium/scripts/panc_xenium_cluster_name.py` (the per-section worker),
`analysis/xenium/scripts/panc_xenium_release.py` (the cohort driver and release writer),
`analysis/xenium/scripts/panc_xenium_external_check.py` (the external check).

**Where the output is.** `/data/xenium_typing/panc_2026-09-17/` — release, per-section work
artifacts and the two 6539 diagnostic trees. `/home` is at 97% with ~110 GB free, so nothing in
this lane writes under `/home`. `data/development/xenium_typing/release_pln_immune_v3/` was not
touched. `data/development/xenium_typing/release/` (the existing 2-donor scANVI release) was read
and never written.

**Latch.** This is a development release. Nothing here is a production typing and no latch moves.

---

## 1. The defect, as it presented

The smoke test at `data/development/xenium_typing/release_panc_v1/_work/6539.calls.parquet`
returned, over 274,531 called cells:

| lineage | cells |
|---|---|
| Exocrine | 176,673 |
| Stromal_Fibroblast_Stellate | 20,856 |
| Vascular | 20,546 |
| Stromal_Pericyte_SMC | 15,707 |
| **Granulocyte** | **15,476** |
| Myeloid | 10,181 |
| Immune_B_Plasma | 6,412 |
| Endocrine | 5,168 |
| Neural | 3,512 |
| **Immune_T** | **0** |
| **Immune_NK_ILC** | **0** |
| (Unavailable_no_nucleus) | 41,755 |
| (Uncalled_below_cluster_qc) | 1,324 |

Zero T cells and zero NK against 15,476 granulocytes.

## 2. Diagnosis — it is two defects, not one, and both were measured

The smoke test was produced by an older worker than the one in the checkout: its report records
`leiden_resolution 1.0`, 19 clusters, no refinement pass and no `score <= 0` guard. Disentangling
the two required running the checked-out worker unchanged, which had never been run.

### 2a. Defect A — the mean over a whole core is biased by core size

The cores are not the same size, and what differs in size is the `cellxgene` supporting genes,
which are broadly expressed. Measured against the 5,101-gene panel:

| core | genes on panel | of which may NAME (`canonical`/`both`) | supporting only |
|---|---|---|---|
| 03_immune_T_cell | 101 | 20 | 81 |
| 06_immune_myeloid | 90 | 31 | 59 |
| 11_endothelial_vascular | 87 | 23 | 64 |
| 04_immune_B_plasma | 58 | 26 | 32 |
| 13_fibroblast_stellate | 39 | 10 | 29 |
| 05_immune_NK_ILC | 31 | 22 | 9 |
| 07_immune_granulocyte | 21 | 19 | 2 |
| 01_pancreas_endocrine | 18 | 16 | 2 |
| 15_neural_glial | 17 | 17 | 0 |
| 12_pericyte_smooth_muscle | 9 | 8 | 1 |
| 02_pancreas_exocrine | 8 | 8 | 0 |

A mean over the whole core is pulled towards the section background in proportion to the
supporting-gene count. Measured on 6539's own 19-cluster partition, at each core's best cluster,
**mean z over naming genes minus mean z over the whole core**:

| core | supporting genes | dilution |
|---|---|---|
| Immune_T | 81 | **+1.033** |
| Vascular | 64 | +0.532 |
| Myeloid | 59 | −0.260 |
| Stromal_Fibroblast_Stellate | 29 | +0.688 |
| Immune_B_Plasma | 32 | +0.344 |
| Immune_NK_ILC | 9 | +0.233 |
| Endocrine | 2 | +0.187 |
| Granulocyte | 2 | +0.086 |
| Stromal_Pericyte_SMC | 1 | +0.190 |

Spearman(supporting-gene count, dilution) = **0.527** over all nine, **0.85** over the eight
excluding Myeloid — whose supporting genes are, as it happens, genuinely myeloid-enriched, which
is why it is the one core the dilution does not touch.

**The consequence was a wrong call.** 6539's cluster `c16` (6,412 cells) is the lymphocyte
cluster. Whole-core scoring: Immune_B_Plasma 2.875, Immune_T 2.628, Immune_NK_ILC 2.531 — B wins
by 0.247 and the whole cluster went out as B. Naming-only scoring on the same numbers: Immune_T
3.661, Immune_B_Plasma 3.219 — T wins by 0.442. **All twenty T naming genes stand out in that
cluster**, with CD3E detected in 37.2% of its cells against 1.1% of the rest of the section, CD8A
30.5% against 0.3%, CD247 14.9% against 0.17%. The cluster does contain B cells too (23 of 26 B
naming genes also stand out); it is a mixed lymphocyte cluster at that resolution, which is
defect B.

### 2b. Defect B — 15,169 of the 15,476 granulocytes rest on one gene in seven cells

6539's cluster `c03`, 15,169 cells, **median depth 64 transcripts per cell**, was called
Granulocyte on a **negative** whole-core score of −0.104 — it won the argmax by being the
least-negative of eleven negatives. The evidence behind that label, from the section's own
gene statistics:

| gene | z | det in c03 | det elsewhere | stands out |
|---|---|---|---|---|
| CCR3 | +2.796 | 0.0005 | 0.0002 | yes |
| CXCR1 | +0.936 | 0.0004 | 0.0002 | no |
| MPO | −0.906 | 0.0001 | 0.0009 | no |
| ELANE | −0.459 | 0.0004 | 0.0027 | no |
| FCGR3B | −0.331 | 0.0004 | 0.0025 | no |
| CSF3R | −0.229 | 0.0010 | 0.0052 | no |
| PRTN3 | −0.101 | 0.0003 | 0.0002 | no |

Every canonical neutrophil gene is **depleted** in that cluster. The entire label rests on CCR3
in about seven of 15,169 cells, below the cluster's own background detection rate of 0.00165.
This is the pLN reviewer's finding — *"a core with few genes can win an argmax on near-zero
evidence"* (`docs/PLN_IMMUNE_XENIUM_2026-09-17.md`) — arriving on the pancreas from the other
side: there the gDNA hazard makes a **large** core easy to anchor by a count rule
(`P(Binomial(n,p) ≥ 2)` grows with `n`, `docs/PLN_IMMUNE_TYPING_HAZARD_2026-09-17.md`); here a
mean rule makes a **small** core easy to win. One underlying defect — neither rule is calibrated
for how many genes a core has — with two faces.

Two guards already in the checked-out worker and never run address this face: the winning score
must be above the partition mean at all, and the standout naming gene must clear the cluster's
**own** background. `c03` fails both.

### 2c. Defect C — one lymphocyte cluster at `resolution 1.0`

At `resolution 1.0` with 19 clusters there is no T cluster in 6539 at all: T, B and NK are one
group, so even a correct score can only return one label for all of them. The checked-out worker
already partitions at `resolution 2.0` and refines every parent of ≥ 2,000 cells; that had never
been run either.

## 3. The fix that was chosen, and the one that was rejected on measurement

**Chosen: score each core on the genes that may NAME its lineage** — `core/panels`' own
`evidence` column, values `canonical` or `both`. The supporting `cellxgene` genes stay in the
embedding, stay in the per-cluster gene statistics and are still written per cluster as
`full_core_score__*`; they no longer vote in the argmax. This removes the size-dependent
dilution because what differed in size was exactly those genes, and it does not invert the bias:
the naming counts run 8 to 31, a far narrower range than the 8 to 101 of the whole cores, and the
score is still a mean, so a core is not rewarded for carrying more naming genes.

**Rejected: normalising each core against its own distribution across the section's clusters.**
Measured on 6539 with a robust centre and MAD across the 19 clusters, it is worse. A core absent
from almost every cluster of the section has a tiny denominator, so its normalised score explodes
(Immune_T 58.3, Vascular 47.8, Endocrine 31.1 on the same clusters where the raw scores are 3.7,
2.6 and 3.6) and the winner becomes whichever core happens to have the **tightest across-cluster
spread** — a property of that section's composition, not evidence. Concretely it moved 6539's
`c04`, a 307-cell cluster carrying six standout granulocyte naming genes (CTSG, ENPP3, GATA2,
IL5RA, PRTN3, RNASE3), to Immune_T, and moved a 767-cell exocrine cluster to B_plasma. It
replaces a bias towards small cores with a bias towards narrow ones. Not adopted.

Nothing else in the method changed. Clustering is still within donor and per section with no
joint embedding; naming is still by enrichment against the other clusters of the same section;
admission is still one standout naming gene of the winning core above the cluster's own
background, and a cluster that fails is uncalled and never reassigned; negative evidence still
redirects a CD3-negative, NK-positive T argmax to NK; detection is still `lognorm >= log1p(1)`
with no raw-count floor; a cell with no nucleus is still excluded from every calculation.

## 4. The effect on 6539, before and after, per lineage

All three columns are the same 317,610 cells of the same section. Column 2 is the checked-out
worker run unchanged (resolution 2.0, refinement on, whole-core score); column 3 is the same run
with naming-only scoring. Isolating them this way is why column 2 exists at all.

| lineage | smoke test (res 1.0, whole-core) | before (res 2.0 + refine, whole-core) | after (res 2.0 + refine, naming-only) |
|---|---|---|---|
| Exocrine | 176,673 | 183,331 | 181,173 |
| Unavailable_no_nucleus | 41,755 | 41,755 | 41,755 |
| Vascular | 20,546 | 20,640 | 20,640 |
| Stromal_Fibroblast_Stellate | 20,856 | 17,805 | 19,026 |
| Stromal_Pericyte_SMC | 15,707 | 16,013 | 16,826 |
| Myeloid | 10,181 | 9,085 | 9,733 |
| Uncalled_no_enrichment | 0 | 9,185 | 5,860 |
| Endocrine | 5,168 | 5,215 | 5,215 |
| **Immune_T** | **0** | **3,519** | **4,229** |
| Neural | 3,512 | 3,716 | 3,631 |
| Uncalled_no_naming_evidence | 0 | 918 | 2,718 |
| **Granulocyte** | **15,476** | **1,818** | **2,118** |
| **Immune_NK_ILC** | **0** | **1,037** | **1,823** |
| Immune_B_Plasma | 6,412 | 2,249 | 1,539 |
| Uncalled_below_cluster_qc | 1,324 | 1,324 | 1,324 |

Read it in two steps, because the two defects do different work:

* **Resolution and the two guards (column 1 → 2)** do most of it. Granulocyte falls from 15,476 to
  1,818 as `c03` loses its negative-score label; Immune_T appears at 3,519 and NK at 1,037 as the
  lymphocyte parent splits into 10 subclusters; 9,185 cells go to `Uncalled_no_enrichment`, almost
  all of them in low-depth parents (`c12`, median depth 47–138; `c14`, median depth 19–39).
* **Naming-only scoring (column 2 → 3)** changes 14 subclusters of 238, 7,770 cells. The largest
  moves: `c20_s07` (710 cells, depth 306) B_Plasma → Immune_T on 18 standout T naming genes
  including CD3E, CD3G, CD4, CD8A; `c17_s00` (2,034 cells, depth 450) Exocrine →
  Stromal_Fibroblast_Stellate — the 39-gene stromal core had been diluted below the 8-gene
  exocrine core; `c22_s06` (85 cells) Neural → Myeloid on eleven myeloid naming genes. Immune_T
  +710, NK +786, Granulocyte +300, B_Plasma −710.

**What the fix does not fix, stated plainly.** Some of the NK and granulocyte gain in column 3 is
in the low-depth parent `c14` (subclusters of median depth 24–32) where the naming-only score
lifts a previously negative score just above zero: `c14_s10` (337 cells, depth 24) → NK on
KIR2DL1/KIT/KLRF1, `c14_s02` (300 cells, depth 28) → Granulocyte on PRG2 alone, `c14_s11` (106
cells, depth 24.5) → Myeloid on LILRA4 alone. These are admissible under the method as stated —
one standout naming gene above the cluster's own background — and they are thin. No threshold was
invented to suppress them. Instead every called cell now carries
`n_naming_standouts_behind_cluster_call` and `cluster_median_depth`, and the per-section report
counts `called_cells_on_one_naming_gene_in_low_depth_cluster`, so a thin call is visible in the
release rather than argued about here.

`c20_s07` is genuinely mixed: 18 T naming genes and 17 B naming genes both stand out in it, and
the margin is 0.367. The cluster frames carry `purity_cell_argmax_agrees`,
`cell_argmax_entropy` and `n_final_labels_in_this_cluster` for exactly this.

---

## 5. The release

`/data/xenium_typing/panc_2026-09-17/release/`

```
release.json                                schema_version 1, release_id
                                            xenium-panc-cluster-name-20260917, roi panc,
                                            21 donors, sha256 per input (21 h5ad + 12 identity
                                            core CSVs + cohort.py)
label_space.json                            11 labels, 4 uncalled states, what each label does
                                            and does not support
calls/donor_id=<d>/roi=panc/data.parquet    24 columns keyed on cell_id, verbatim
per_donor_summary.csv
lineage_counts.csv
coherence_genes.csv
coherence_where_the_evidence_went.csv
external_check_vs_scanvi.json
```

`donor_id` and `roi` are the hive partition keys and are **not** repeated as columns, so `calls/`
reads as one dataset. **Read `donor_id` as a string** — every donor in this cohort is all digits,
so an inferring reader turns it into an integer and `112` stops matching the donor it names. The
exact recipe is in `release.json` under `read_calls_with`:

```python
import pyarrow as pa, pyarrow.dataset as ds
part = ds.partitioning(pa.schema([("donor_id", pa.string()), ("roi", pa.string())]), flavor="hive")
calls = ds.dataset("/data/xenium_typing/panc_2026-09-17/release/calls",
                   format="parquet", partitioning=part)   # 17,490,305 rows
```

`nucleus_count` is a **nullable** Int32: where the nucleus rule was unmeasurable it is `NA` for
every cell of that donor, never `0` and never a negative sentinel, and `nucleus_measured` is
`False` beside it.

Per-section work artifacts (cluster frames, subcluster frames, gene statistics, reports) are in
`/data/xenium_typing/panc_2026-09-17/_work/`; the two 6539 diagnostic trees are `_diag_before/`
and `_diag_after/`. Nothing under `/home` was written and
`data/development/xenium_typing/release/` and `.../release_pln_immune_v3/` were left alone.

**21 sections, 21 donors, 17,490,305 cells, 15,318,868 called (87.6%), 775,234 excluded as having
no nucleus.** 6457 is excluded via `phenocycler.cohort.DONOR_EXCLUSIONS`; 6579 is in, reinstated
2026-09-17. The h5ads carry one pancreas section per donor, so section = donor here.

### Per donor

| donor | sex | status | batch | cells | median tx/cell | clusters / subclusters | kNN holdout | called | called of nucleated |
|---|---|---|---|---|---|---|---|---|---|
| 6347 | M | AAb+ | 3 | 1,393,487 | 161 | 25 / 331 | 0.889 | 1,230,053 | 0.883 |
| 6374 | F | ND | 6 | 1,012,696 | 420 | 22 / 210 | 0.908 | 868,022 | 0.917 |
| 6380 | F | T1D | 5 | 689,516 | 221 | 23 / 220 | 0.885 | 606,873 | 0.896 |
| 6436 | M | T1D | 10 | 692,251 | 420 | 17 / 198 | 0.906 | 613,481 | 0.922 |
| 6442 | M | ND | 11 | 956,008 | 377 | 26 / 232 | 0.889 | 840,811 | 0.930 |
| 6461 | M | ND | 7 | 812,242 | 438 | 26 / 232 | 0.883 | 660,181 | 0.902 |
| 6476 | F | ND | 4 | 1,109,281 | 493 | 25 / 245 | 0.899 | 924,171 | 0.882 |
| 6482 | F | ND | 9 | 810,102 | 474 | 25 / 242 | 0.896 | 710,510 | 0.920 |
| 6493 | M | ND | 8 | 655,839 | 491 | 24 / 244 | 0.899 | 542,573 | 0.904 |
| 6521 | M | AAb+ | 8 | 495,617 | 518 | 27 / 237 | 0.890 | 397,166 | 0.861 |
| 6523 | F | T1D | 5 | 687,186 | 296 | 21 / 228 | 0.883 | 580,219 | 0.907 |
| 6533 | F | T1D | 2 | 1,048,928 | 421 | 32 / 328 | 0.911 | 990,951 | 0.964 |
| 6539 | M | ND | 10 | 317,610 | 330 | 30 / 238 | 0.907 | 265,953 | 0.964 |
| 6551 | M | T1D | 9 | 760,684 | 620 | 32 / 346 | 0.919 | 690,534 | 0.975 |
| 6566 | M | T1D | 7 | 831,866 | 707 | 32 / 317 | 0.920 | 728,005 | 0.941 |
| 6579 | F | T1D | 6 | 665,135 | 409 | 23 / 243 | 0.922 | 619,490 | 0.931 |
| 6580 | M | T1D | 3 | 931,227 | 509 | 25 / 260 | 0.915 | 846,117 | 0.909 |
| 6591 | M | T1D | 1 | 1,080,047 | 410 | 24 / 255 | 0.904 | 949,091 | 0.902 |
| 6593 | M | T1D | 11 | 1,059,148 | 540 | 31 / 314 | 0.910 | 958,867 | 0.928 |
| 112 | F | ND | 2 | 546,910 | 434 | 26 / 243 | 0.899 | 473,207 | 0.912 |
| 115 | F | ND | 1 | 934,525 | 430 | 24 / 246 | 0.900 | 822,593 | 0.927 |

`kNN holdout` is the measured cost of learning the Leiden partition on a 60,000-cell subsample and
carrying it outward by a 15-nearest-neighbour vote: the fraction of held-out subsample cells the
classifier puts back in their own Leiden cluster. It is reported, not assumed.

### Per donor, per lineage, including uncalled

| donor | Endocrine | Exocrine | Immune_T | NK_ILC | B_Plasma | Myeloid | Granulocyte | Vascular | Pericyte_SMC | Fibro_Stellate | Neural | Unc: below QC | Unc: no enrichment | Unc: no naming | Unavail: no nucleus |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 6347 | 43,427 | 727,943 | 38,979 | 1,702 | 45,204 | 157,303 | 11,002 | 69,556 | 17,799 | 110,403 | 6,735 | 23,271 | 108,466 | 31,697 | unmeasured |
| 6374 | 22,757 | 632,190 | 3,125 | 1,488 | 28,900 | 32,857 | 16,386 | 50,837 | 14,157 | 62,513 | 2,812 | 2,433 | 69,255 | 7,316 | 65,670 |
| 6380 | 9,850 | 403,467 | 2,706 | 2,046 | 449 | 58,906 | 12,860 | 33,497 | 10,976 | 65,022 | 7,094 | 5,083 | 58,368 | 7,098 | 12,094 |
| 6436 | 21,512 | 422,254 | 4,153 | 103 | 10,323 | 23,275 | 28,639 | 36,433 | 14,404 | 50,035 | 2,350 | 4,555 | 38,066 | 9,341 | 26,808 |
| 6442 | 10,818 | 632,571 | 11,327 | 433 | 2,490 | 76,029 | 10,846 | 32,008 | 7,079 | 56,575 | 635 | 3,207 | 55,133 | 5,471 | 51,386 |
| 6461 | 14,313 | 481,755 | 1,565 | 269 | 25,906 | 21,337 | 4,694 | 38,687 | 7,667 | 57,887 | 6,101 | 1,199 | 57,183 | 13,367 | 80,312 |
| 6476 | 30,547 | 645,566 | 3,662 | 461 | 10,959 | 77,471 | 3,656 | 50,706 | 13,562 | 84,251 | 3,330 | 2,080 | 117,595 | 3,773 | 61,662 |
| 6482 | 20,361 | 529,293 | 2,133 | 865 | 8,696 | 23,076 | 4,965 | 44,554 | 15,108 | 56,266 | 5,193 | 1,289 | 46,060 | 14,097 | 38,146 |
| 6493 | 13,417 | 385,601 | 1,674 | 550 | 24 | 41,967 | 3,521 | 33,195 | 9,289 | 48,032 | 5,303 | 1,447 | 46,591 | 9,757 | 55,471 |
| 6521 | 9,848 | 269,697 | 3,695 | 135 | 16,306 | 17,636 | 2,247 | 37,204 | 8,193 | 29,833 | 2,372 | 876 | 52,214 | 11,104 | 34,257 |
| 6523 | 15,390 | 338,386 | 8,025 | 902 | 76,551 | 22,996 | 5,273 | 43,919 | 13,077 | 51,363 | 4,337 | 3,464 | 45,032 | 11,236 | 47,235 |
| 6533 | 28,359 | 590,397 | 53,560 | 2,198 | 30,294 | 64,291 | 11,458 | 80,250 | 27,496 | 95,672 | 6,976 | 2,641 | 12,789 | 21,703 | 20,844 |
| 6539 | 5,215 | 181,173 | 4,229 | 1,823 | 1,539 | 9,733 | 2,118 | 20,640 | 16,826 | 19,026 | 3,631 | 1,324 | 5,860 | 2,718 | 41,755 |
| 6551 | 12,754 | 497,145 | 23,303 | 457 | 7,323 | 23,011 | 2,912 | 49,314 | 22,504 | 47,236 | 4,575 | 714 | 12,883 | 4,071 | 52,482 |
| 6566 | 12,678 | 463,592 | 11,875 | 2,426 | 8,245 | 83,092 | 37,700 | 35,980 | 13,523 | 55,766 | 3,128 | 2,349 | 30,334 | 13,354 | 57,824 |
| 6579 | 25,535 | 357,025 | 30,436 | 5,203 | 16,353 | 49,158 | 35,729 | 39,886 | 8,156 | 48,766 | 3,243 | 4,114 | 38,645 | 2,886 | unmeasured |
| 6580 | 31,385 | 606,561 | 6,038 | 529 | 1,566 | 33,080 | 4,099 | 73,813 | 17,243 | 63,864 | 7,939 | 2,451 | 57,733 | 24,926 | unmeasured |
| 6591 | 47,506 | 647,444 | 12,188 | 700 | 376 | 25,855 | 4,556 | 83,671 | 26,720 | 91,728 | 8,347 | 1,902 | 100,936 | 398 | 27,720 |
| 6593 | 17,691 | 597,958 | 15,097 | 1,762 | 855 | 99,159 | 49,019 | 62,018 | 20,073 | 92,428 | 2,807 | 1,584 | 60,864 | 11,522 | 26,311 |
| 112 | 20,241 | 328,194 | 1,759 | 579 | 38 | 30,048 | 3,154 | 27,395 | 8,498 | 52,589 | 712 | 1,169 | 28,823 | 15,940 | 27,771 |
| 115 | 49,575 | 549,248 | 782 | 1,601 | 229 | 65,885 | 6,598 | 57,443 | 17,113 | 69,048 | 5,071 | 2,854 | 43,481 | 18,111 | 47,486 |

`unmeasured` in the last column is exactly that. The h5ads carry no nucleus geometry, so the only
handle on the nucleus rule is the vendor bundle's `cells.parquet`, matched by donor **and** by
subset containment of the h5ad's own cell ids — never by filename order. For **6347, 6579 and
6580** that bundle is not on this machine. The batch trees themselves are present; only their
**pancreas** output directories are missing, while their pLN siblings are there:

| donor | missing directory (under `/data/xenium/HuPrime5kPanTissue/6KTRAF/`) | siblings that are present |
|---|---|---|
| 6580 | `MM037_Xenium_nPODBatch3/output-XETG00298__0041139__Panc__20241025__204655` | `…0041139__pLN…`, `…0041140__pLN…` |
| 6347 | `MM037_Xenium_nPODBatch3/output-XETG00298__0041140__Panc_pLN__20241025__204655` | as above |
| 6579 | `MM040_Xenium_nPODBatch6/output-XETG00298__0041035__Pan__20241219__201934` | `…0041035__pLN1/2/3…`, `…0041032__Pan…` |

So no local `cells.parquet` is a superset of those h5ads' cell ids and the nucleus rule could not
be applied. It is recorded as unmeasurable in each of those
reports (`nucleus_rule_measurable: false`) — not as zero, and not as passed. Their counts above
therefore still contain cells that may have no nucleus, and that is why
`called_fraction_of_all_cells` equals `called_fraction_of_nucleated_cells` for exactly those three.

**This needs the maintainer.** A cell with no nucleus is excluded from every calculation, and for
three of 21 donors that rule is currently unenforceable here. Either those three pancreas output
directories are fetched, or `xenium_paths.csv` is corrected to the bundle those h5ads were actually
built from. Until then, 6347, 6579 and 6580 should not be pooled with the other eighteen for any
per-cell rate. In the other eighteen the excluded fraction runs 1.8%–13.2% of all cells, so the gap
is not negligible.

Cohort totals: Exocrine 10,287,460 (58.8%) · Stromal_Fibroblast_Stellate 1,308,303 (7.5%) ·
Myeloid 1,036,165 (5.9%) · Vascular 1,001,006 (5.7%) · Endocrine 463,179 (2.6%) ·
Stromal_Pericyte_SMC 309,463 (1.8%) · Immune_B_Plasma 292,626 (1.7%) · Granulocyte 261,432 (1.5%) ·
Immune_T 240,311 (1.4%) · Neural 92,691 (0.5%) · Immune_NK_ILC 26,232 (0.15%); uncalled:
no_enrichment 1,086,311 (6.2%), no_naming_evidence 239,886 (1.4%), below_cluster_qc 70,006 (0.4%);
no nucleus 775,234 (4.4%).

**Every one of the 21 sections has non-zero Immune_T and non-zero Immune_NK_ILC.**

### What the Endocrine and Exocrine labels do and do not support

INS, GCG, SST and PPY are absent from the 5,101-gene panel — measured per section and recorded in
every report, not assumed. Endocrine resolves only through the curated surrogate panel
`01_pancreas_endocrine` (18 genes, 16 of them naming: ABCC8, ADCYAP1, ARX, FEV, INSM1, ISL1,
KCNJ11, MAFB, NEUROD1, PAX6, PCSK1, PDX1, SCG2, SCG5, SLC30A8 and others). The call supports *this
cluster is endocrine*. It supports **no alpha, beta or delta statement**, and none is made anywhere
in this release.

Exocrine is **one label**. PRSS1, CPA1, CELA3A, CTRB1, CTRB2, AMY2A and REG1A are all absent from
the panel, measured per section. **No acinar and no ductal claim is made or supported**, and
`16_epithelial_general` is an embedding extra that is measured, reported per cluster as
`extra_score__epithelial_general`, and never a lineage — KRT and EPCAM are shared by ductal,
acinar and endocrine cells.

## 6. Controls

### The Y-probe control, on female donors only

The existing 6539 report recorded `y_probe_detect_frac_any 0.260` and 6539 is **male**, so that
number was a measurement of 6539's Y expression, not a control. Sex was resolved per donor from
each h5ad's own `obs['sex']`: nine female donors (6374, 6380, 6476, 6482, 6523, 6533, 6579, 112,
115) and twelve male.

The floor group is TBL1Y / USP9Y / EIF1AY — zero-template in a female donor and, unlike DDX3Y and
KDM5D, not cross-hybridising with a highly expressed X paralog. Every detection of those in a
female donor is false, so the rate reads the background this detection predicate admits.

| | cell-weighted detect frac (any floor probe) |
|---|---|
| **female donors (the control)** | **0.000196 — 0.0196%** |
| target | ≤ ~0.3% |
| male donors (not a control) | 0.1479 |

Per female section: 6374 0.00039 · 6380 0.00010 · 6476 0.00030 · 6482 0.00016 · 6523 0.00009 ·
6533 0.00008 · 6579 0.00019 · 112 0.00008 · 115 0.00025. Every one is an order of magnitude inside
the target. The male rate is what a real Y template looks like under the same predicate and is
reported for contrast, never used as a background estimate.

### Depth

Depth, not gDNA, is the constraint on a 5K-panel call for most of this cohort.

> **CORRECTED 2026-09-18.** This sentence used to end "and no section here is described as
> gDNA-hazardous". **Two sections in this release are exactly the flagged pair.**
> `docs/XENIUM_GDNA_TRIAGE_2026-09-17.md` measures `0080870` panc (donor **6533**) at 0.26100,
> **57.5x** the tuned rate, and `0081366` panc (donor **112**) at 0.30600, **67.4x**;
> `docs/XENIUM_PIPELINE.md` §0 puts them at roughly **76x** the false positives per cell of the
> sample every constant in this lane was tuned on. Both donors are typed in this release (see
> the per-donor table above) with **no gDNA correction applied and no flag on the rows**. Any
> call on 6533 or 112 from this release has to be read with that in front of it, and a
> cohort-pooled number that includes them carries it silently. Median transcripts per cell, lowest first: 6347 161 · 6380 221 · 6523 296 ·
6539 330 · 6442 377 · 6579 409 · 6591 410 · 6374 420 · 6436 420 · 6533 421 · 115 430 · 112 434 ·
6461 438 · 6482 474 · 6493 491 · 6476 493 · 6580 509 · 6521 518 · 6593 540 · 6551 620 · 6566 707.

**No section is flagged low depth** at the section level (the flag fires at a median of 100 or
below; the lowest is 6347 at 161). Low depth is nonetheless real *within* sections: it concentrates
in particular clusters, and those are where the thin calls are. Every called cell carries
`cluster_median_depth`, so a call made on a 24-transcript cluster is identifiable in the release.

### How decided each call was

The method admits the argmax on evidence — one standout naming gene of the winning core, above the
cluster's own background — and **no margin threshold operates anywhere**. That means some admitted
calls are decided and some are not, and both go out. Rather than invent a cutoff, the release
carries the decidedness on every called cell: `cluster_score_margin`, `cluster_naming_standouts`,
`cluster_purity_cell_argmax_agrees`, `cluster_cell_argmax_entropy`,
`n_naming_standouts_behind_cluster_call`, `cluster_median_depth`.

Of 15,318,868 called cells: **440,167 (2.9%) were called by a cluster that beat its runner-up by
less than 0.05**, and **51,434 (0.34%) by a cluster that rested on a single naming gene and had a
median depth below 50 transcripts**. The largest single example is `6436_c04_s07`: 10,241 cells,
median depth 368, called Immune_B_Plasma on FKBP11 and XBP1 with a margin of 0.008 and 0.2% of its
own cells agreeing by their own core argmax. It is in the release with all six of those numbers
attached.

> **CORRECTED 2026-09-18.** This passage used to read: "`n_cd3_negative_nk_redirects` is **0**
> across all 21 sections: no cluster this time had an Immune_T argmax with neither CD3E nor CD3G
> standing out while NK genes did. The rule is live and simply did not fire." **The rule fired.**
> The published 0 is a counter reading the wrong partition, not a measurement of the rule.
>
> `analysis/xenium/scripts/panc_xenium_cluster_name.py:885` computes
> `n_clusters_cd3_negative_nk_redirect` as `int(sum(main["redirected"]))` — the **parent-cluster**
> decision and only that — and `panc_xenium_release.py:184` copies it into `release.json` under the
> plural name `n_cd3_negative_nk_redirects`. The redirect is evaluated again at the **sub-cluster**
> level (`fdec["redirected"]`, written onto every sub-cluster row at line 821), and the sub-cluster
> partition is where it fired. Nothing in `release.json` counts that partition, so the release
> reports 0 while carrying the redirected cells.

`n_cd3_negative_nk_redirects` is **0** across all 21 sections **because it counts parent clusters
only**, and at the parent level it is a true zero: 0 of 540 parent clusters across the 21 sections
carry `cd3_negative_nk_redirect = True`. At the sub-cluster level the rule fired on **4 of 5,409
sub-clusters**. Measured 2026-09-18 over every section's
`/data/xenium_typing/panc_2026-09-17/_work/*.clusters.parquet` and `*.subclusters.parquet`:

| sub-cluster | donor | cells | parent lineage | final lineage | standout naming genes | margin | median depth |
|---|---|---|---|---|---|---|---|
| `6482_c18_s07` | 6482 | **215** | Exocrine | **Immune_NK_ILC** (`T_NK`) | `EOMES;GATA3;KLRC1;PRF1` | 0.075 | 91 |
| `6347_c20_s03` | 6347 | 211 | Neural | `Uncalled_no_enrichment` | — | 0.126 | 26 |
| `6380_c11_s11` | 6380 | 763 | `Uncalled_no_enrichment` | `Uncalled_no_enrichment` | — | 0.017 | 29 |
| `6523_c19_s04` | 6523 | 109 | Granulocyte | `Uncalled_no_enrichment` | — | 0.280 | 138 |

Three of the four were redirected and then failed the enrichment gate — their winning
`Immune_NK_ILC` score sat below the partition mean — and ended `Uncalled_no_enrichment`, so the
redirect changed no label there. **One produced a call.** `6482_c18_s07` is 215 cells that the
parent cluster held as Exocrine and that the redirect moved to `Immune_NK_ILC`, common lineage
`T_NK`, on four standout NK naming genes with neither CD3E nor CD3G standing out. All 215 are in
the release at `calls/donor_id=6482/roi=panc/data.parquet` carrying
`reason = cd3_negative_nk_positive_redirect`; that reason string appears on **215 cells
release-wide and on no others**.

Read that call with its own decidedness attached, as §6 requires of every call: its
`purity_cell_argmax_agrees` is **0.042** and its `cell_argmax_top_label` is **Exocrine at 0.674** —
the sub-cluster is named NK by cluster-level enrichment while most of its cells still argmax
Exocrine by their own cores. 215 cells is 0.028% of donor 6482's 770,667 clustered-and-called
cells (810,102 in the h5ad, 38,146 excluded for no nucleus, 1,305 below cluster QC). Nothing
here is a verdict on whether those cells are NK; what is corrected is the claim that the rule never
fired.

## 7. Internal coherence — where the naming evidence actually went

**This is a section compared against itself. It is not a check against any prior biological
finding, no expected T-cell fraction exists anywhere in this lane, and nothing was tuned toward
one.**

The coarse test — a lineage receiving zero cells while its own naming genes are detected — passes:
`sections_failing: none`. Section-wide detection of the coherence genes, all on panel in all 21
sections, ranges CD3E 0.0025–0.0517, CD3G 0.0005–0.0062, CD4 0.0117–0.0835, CD8A 0.0014–0.0196,
PTPRC 0.0997–0.2538, NKG7 0.0011–0.0179, KLRD1 0.0003–0.0024, MS4A1 0.0006–0.0329, CD68
0.0160–0.1775 (`coherence_genes.csv`).

That test is too coarse on its own, so a finer one was run
(`analysis/xenium/scripts/panc_xenium_coherence.py`, output
`coherence_where_the_evidence_went.csv`): **of the cells that detect a naming gene, what fraction
sit in a subcluster called that lineage?**

| section | CD3E → Immune_T | MS4A1 → B_Plasma | CD68 → Myeloid | NKG7 → NK_ILC | FCGR3B → Granulocyte |
|---|---|---|---|---|---|
| 6347 | 0.681 | 0.841 | 0.551 | 0.001 | 0.003 |
| 6374 | 0.195 | 0.058 | 0.397 | 0.028 | 0.098 |
| 6380 | 0.294 | 0.170 | 0.575 | 0.002 | 0.107 |
| 6436 | 0.162 | 0.047 | 0.228 | 0.020 | 0.346 |
| 6442 | 0.456 | 0.044 | 0.498 | 0.049 | 0.024 |
| 6461 | 0.213 | 0.063 | 0.371 | 0.138 | 0.028 |
| 6476 | 0.436 | 0.051 | 0.539 | 0.137 | 0.059 |
| 6482 | 0.271 | 0.056 | 0.444 | 0.004 | 0.085 |
| 6493 | 0.248 | 0.029 | 0.475 | 0.163 | 0.012 |
| 6521 | 0.487 | 0.054 | 0.497 | 0.045 | 0.029 |
| 6523 | 0.271 | 0.211 | 0.506 | 0.000 | 0.020 |
| 6533 | 0.691 | 0.782 | 0.343 | 0.173 | 0.028 |
| 6539 | 0.312 | 0.403 | 0.465 | 0.276 | 0.006 |
| 6551 | 0.657 | 0.599 | 0.329 | 0.093 | 0.052 |
| 6566 | 0.272 | 0.514 | 0.348 | 0.121 | 0.197 |
| 6579 | 0.520 | 0.653 | 0.264 | 0.104 | 0.183 |
| 6580 | 0.482 | 0.511 | 0.380 | 0.154 | 0.206 |
| 6591 | 0.694 | 0.319 | 0.181 | 0.086 | 0.016 |
| 6593 | 0.407 | 0.053 | 0.429 | 0.109 | 0.553 |
| 112 | 0.223 | 0.029 | 0.501 | 0.169 | 0.064 |
| 115 | 0.054 | 0.051 | 0.515 | 0.002 | 0.043 |

**This is the honest residual, and it is a resolution limit, not a scoring one.** Traced in the
three weakest sections:

* **115** — CD3E detected in 6,001 cells; 326 of them sit in an Immune_T subcluster, 4,173 in
  Exocrine subclusters. The T subclusters that do exist are unambiguous (`HDL115_c00_s07`, 491
  cells, CD3E in 55.2% of them, z = 14.8, nineteen T naming genes above background). The CD3E
  outside them sits in large, deep exocrine subclusters at 1.5–3.8% detection
  (`HDL115_c06_s00`, 17,520 cells, median depth 822, CD3E 1.7%, z = 0.10, not standing out).
* **6436** — CD3E in 9,042 cells, 1,463 in an Immune_T subcluster. `6436_c01_s00`: 2,403 cells,
  median depth 722, **CD3E detected in 28.2% and standing out at z = 2.75**, yet called Exocrine
  on three exocrine naming genes with a margin of 0.38. That subcluster contains T cells and did
  not separate from exocrine.
* **6493** — CD3E in 2,659 cells, 660 in an Immune_T subcluster, and a further 291 in
  `6493_c22_s05` (550 cells, CD3E 52.9%), which was called Immune_NK_ILC over Immune_T by 0.25.

The method says cells inherit the cluster call and per-cluster heterogeneity is reported so a mixed
cluster is visible. That is exactly what is happening: where T cells are sparse and embedded in
exocrine tissue they do not form their own subcluster at this partition, and the cells inherit the
majority label. Nothing was changed to move those numbers. `cell_core_argmax` and
`cell_agrees_with_cluster_call` are on every cell so this is measurable downstream rather than
invisible; **the next lever on this lane is resolution inside the large exocrine parents, not the
score.**

## 8. External check — the new calls against the existing 2-donor scANVI release

Nothing in this lane was fitted to, tuned against or gated on
`data/development/xenium_typing/release/`. It was read once, after the release was written.

The two lanes type different cells by different means from different h5ads of the same sections,
and their label spaces differ: the old space has no granulocyte and no NK label and collapses both
stromal cores into `Stromal`. The new space was mapped onto the old one (NK → Immune_T,
Granulocyte → Myeloid, both stromal cores → Stromal), which is the only direction that loses no
distinction the old release ever made.

| | 6436 | 6442 |
|---|---|---|
| cell ids in both | 692,251 | 956,008 |
| comparable (called here, not `Unavailable` there) | 581,148 | 810,726 |
| **agreement** | **0.814** | **0.832** |
| Endocrine | 0.964 | 0.904 |
| Exocrine | 0.841 | 0.932 |
| Vascular | 0.918 | 0.287 |
| Myeloid | 0.827 | 0.831 |
| Immune_T | 0.646 | 0.575 |
| Stromal | 0.669 | 0.148 |
| Neural | 0.422 | 0.477 |
| Immune_B | 0.006 | 0.068 |

The largest single disagreement is T cells in 6436: the old lane calls 42,603 cells Immune_T where
this one calls 4,153, and 29,807 cells are Exocrine here and Immune_T there. **Neither number is a
verification target and this comparison does not adjudicate between the lanes.** What can be
measured is how much direct T evidence each set of calls carries, on the same cells, at raw
count ≥ 1:

| set of cells | n | CD3E | CD3G | CD2 | CD8A | any of the six | PTPRC |
|---|---|---|---|---|---|---|---|
| 6436, called T **here** | 4,153 | 0.352 | 0.041 | 0.102 | 0.367 | **0.621** | 0.450 |
| 6436, called T **there** | 42,603 | 0.069 | 0.007 | 0.012 | 0.043 | **0.120** | 0.195 |
| 6436, Exocrine here & T there | 29,807 | 0.051 | — | — | — | **0.081** | 0.177 |
| 6436, whole section | 692,251 | 0.013 | — | — | — | **0.023** | 0.180 |
| 6442, called T **here** | 11,327 | 0.477 | 0.116 | 0.236 | 0.436 | **0.813** | 0.696 |
| 6442, called T **there** | 12,038 | 0.287 | 0.078 | 0.154 | 0.283 | **0.506** | 0.492 |
| 6442, whole section | 956,008 | 0.013 | — | — | — | **0.026** | 0.219 |

Stated as a measurement and nothing more: on both donors the cells called T by this lane detect T
evidence several times more often than the cells called T by the other lane, and the 29,807 cells
the two lanes disagree about in 6436 detect it at 0.081 against a whole-section 0.023 — above
background, well below either lane's T calls. This says the disagreement is about cells with weak
and mixed evidence. It does not say which lane is right about them, and no call in this release was
moved because of it.

## 9. What is not claimed

* No acinar or ductal call. No alpha, beta or delta call. Exocrine is one label.
* No spatial, per-islet or disease-progression statement of any kind is made here; this document is
  the typing and its defects.
* No result in this document is a verification target. In particular, the T-cell fraction is
  reported per donor and is not tuned, gated or compared against any expected value.
* The latch is untouched. This is a development release.

## 10. How to reproduce it

```bash
cd /home/smith6jt/IO60panc2nd/Islet-Explorer-Senior
.venv/bin/python analysis/xenium/scripts/panc_xenium_release.py \
    --work /data/xenium_typing/panc_2026-09-17/_work \
    --release /data/xenium_typing/panc_2026-09-17/release \
    --score-on naming --jobs 4
.venv/bin/python analysis/xenium/scripts/panc_xenium_coherence.py
.venv/bin/python analysis/xenium/scripts/panc_xenium_external_check.py --evidence
```

One section alone, and the whole-core score for the before/after:

```bash
.venv/bin/python analysis/xenium/scripts/panc_xenium_cluster_name.py \
    --section 6539 --work /data/xenium_typing/panc_2026-09-17/_diag_after --score-on naming
.venv/bin/python analysis/xenium/scripts/panc_xenium_cluster_name.py \
    --section 6539 --work /data/xenium_typing/panc_2026-09-17/_diag_before --score-on all
```

21 sections on 4 processes took 26 minutes wall; the slowest section was 6347 at 381 s (1,393,487 cells).
