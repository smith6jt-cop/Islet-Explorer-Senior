# pLN Xenium — immune typing, and the T/NK split made by the gene panel

**The maintainer's ruling, 2026-09-17.** *"You say NK/ILC cannot be named when all you have to do is
separate the cells based on CD8 and then apply the gene panel to those cells."*

They are right and the earlier reasoning was wrong. It argued that NK could not be named because no
single probe on either panel is NK-specific. **Negative evidence is evidence.** NK/ILC is the
CD3-negative half of a class the panels already carry as one: PhenoCycler's `lineages.py` holds
`T_NK: (CD3e, CD8, CD4, FOXP3, CD56)`, and `core/panels/data/identity_core/05_immune_NK_ILC.csv`
holds 31 genes, 22 of which may name the lineage and all 22 of which are on the 5K panel. The coarse
class is what the panel carries; the split is made by the presence/absence pattern inside it. That
is what this build does, on both units where the data supports it.

| | |
|---|---|
| Producer | `analysis/xenium/scripts/pln_immune_typing_v3.py` |
| Release | `data/development/xenium_typing/release_pln_immune_v3/` |
| `release_id` | `xenium-pln-immune-tnk-split-dev-20260917`, `schema_version` 1, `roi` `pln_1` |
| Cells | **6,007,789** in 22 sections, **20 donors**; **5,882,198** partitioned |
| T_NK compartment | **1,659,255 cells**, split per cell and per sub-cluster |
| Latch | development; nothing here is a release of record |

`release_pln_immune_v2/` and `release_pln_immune/` are left in place and are superseded, not
overwritten. `data/development/typing/pln/clusters/xenium/` — the shared clustering track — was
**read, never written**.

---

## 1. What the T/NK split is

Two units, two questions, both reported. Neither replaces the other and **neither was tuned to the
other or to the PhenoCycler side.**

### Per cell — the floor

For a cell in the T_NK compartment:

```
CD3E or CD3G detected     ->  T, then:  FOXP3 -> T_Treg
                                        CD8A/CD8B and not CD4 -> T_CD8
                                        CD4 and not CD8A/CD8B -> T_CD4
                                        both -> T_CD4_CD8_both
                                        none of them -> T_subset_unresolved
neither detected, and >= 2 NK-exclusive naming genes detected  ->  NK_ILC
neither detected, fewer than 2                                 ->  T_NK_unresolved
```

Every `NK_ILC` cell here is **individually corroborated**: on that one cell, two NK naming genes
present and both CD3 chains absent. At a median of ~200 transcripts per cell most real NK cells will
not clear that, so this is a **floor**, and the size of the undercount is not something a per-cell
rule can measure.

### Per T_NK sub-cluster — the population estimate

The T_NK compartment of each section is sub-partitioned (same algorithm, same seed derivation as the
coarse partition) and each sub-cluster is named against the **other T_NK sub-clusters of that same
section**:

```
NK_ILC  <=>  CD3E or CD3G STANDS DOWN  (z <= -1 against the compartment AND detected less often
                                        than in the rest of it)
        AND  >= 2 NK-exclusive naming genes STAND UP
        AND  the fraction of its cells detecting at least one of them exceeds that fraction in
             the rest of the compartment
```

This is the unit the PhenoCycler pLN lane also calls at, which is why the cross-modality comparison
in §7 uses it.

> **Two rules were tried and refuted before this one, and both are recorded because the numbers they
> produced are on the record too.**
>
> 1. *"CD3 does not stand out"* as the negative. Inside a compartment that is mostly T, not standing
>    out is the default state of nine sub-clusters in ten. It returned **35.0% NK on donor 6347**.
>    That is a rule failing, not a finding. **Standing down** — the exact mirror of standout, same
>    contrast, same units, opposite direction — is what "CD3-negative" has to mean when every number
>    is relative.
> 2. *Two standout NK genes, with no set-level check.* Asked of two genes out of nineteen, that
>    passes in sub-clusters whose cells carry **less** NK evidence than the compartment average —
>    not thin evidence, contradicted evidence. Nine such sub-clusters (28,361 cells) were admitted
>    and are named in §5. The set-level detection contrast removes them. It is not a new threshold:
>    it is the detection half of the standout rule, asked of the set instead of gene by gene.

### Which genes may name which lineage

Read off `core/panels`' own `evidence` column — `canonical` or `both` may name, `cellxgene`
("marks lineage (CL)") describes only.

- **NK-exclusive naming set, 19 genes, all on panel:** `CD244 EOMES FCGR3A GATA3 IL1RL1 KIR2DL1
  KIR2DL4 KIR3DL1 KIT KLRB1 KLRC1 KLRC2 KLRD1 KLRF1 NCAM1 NCR1 NCR3 RORC TBX21`.
  `KLRG1`, `KLRK1` and `PRF1` are naming genes of the NK core **but also appear in the T core**, so
  they are dropped from the NK evidence: a gene whose presence is equally explained by a T cell
  cannot name NK.
- **T naming set, 20 genes, all on panel:** `CD2 CD247 CD28 CD3E CD3G CD4 CD5 CD6 CD7 CD8A CD8B
  CTLA4 FOXP3 IKZF2 IL2RA ITK LAT LCK THEMIS ZAP70`.
- **Not used for the gate, and why.** `CD3D` is absent from the 5K panel — measured, never
  substituted. **`CD247` (CD3-ζ) is deliberately not used**: NK cells carry CD3-ζ, so it would
  collapse the distinction being drawn. `TRAC` is genuinely T-restricted and is on the panel but is
  not in the panel authority's T core, so it is carried as a per-cell column and a sensitivity, never
  as the gate.
- `GATA3`, `IL1RL1`, `KIT` and `RORC` carry the core's own ILC `subtype_hint`. The identity core is
  one combined NK/ILC core, so **the label is combined**; `n_nk_like_detected` and
  `n_ilc_hint_detected` are written per cell as description, not as a call.

---

## 2. The result, per donor

`NK_ILC (cell)` is the per-cell floor; `NK_ILC (sub-cl.)` is the population estimate. Both are
counts of cells.

| donor | cells | % immune | T_NK | T_CD4 | T_CD8 | T_Treg | CD4+CD8 | T subset unres. | **NK_ILC (cell)** | **NK_ILC (sub-cl.)** | % of T_NK | T_NK unres. |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 112 | 221,578 | 49.1 | 55,328 | 3,796 | 12,716 | 1,445 | 1,982 | 15,383 | **1,534** | **3,600** | 6.51 | 18,472 |
| 115 | 247,399 | 49.7 | 55,906 | 4,436 | 9,992 | 1,064 | 704 | 18,885 | **179** | 0 | 0.00 | 20,646 |
| 6347 | 79,522 | 49.0 | 14,459 | 1,778 | 2,007 | 644 | 428 | 4,384 | **338** | **842** | 5.82 | 4,880 |
| 6374 | 164,364 | 39.6 | 34,251 | 1,779 | 3,976 | 519 | 499 | 10,350 | **636** | 0 | 0.00 | 16,492 |
| 6380 | 264,010 | 56.6 | 79,213 | 3,132 | 7,272 | 2,211 | 449 | 22,571 | **369** | 0 | 0.00 | 43,209 |
| 6436 | 370,852 | 59.7 | 134,562 | 13,153 | 29,755 | 6,869 | 6,704 | 39,386 | **737** | **10,078** | 7.49 | 37,958 |
| 6442 | 176,137 | 66.7 | 69,808 | 17,109 | 15,990 | 1,415 | 7,629 | 15,618 | **2,090** | **3,218** | 4.61 | 9,957 |
| 6461 | 109,476 | 41.0 | 19,481 | 918 | 1,613 | 848 | 106 | 4,590 | **578** | **1,060** | 5.44 | 10,828 |
| 6476 | 290,429 | 59.4 | 53,152 | 3,379 | 7,594 | 700 | 598 | 17,642 | **447** | 0 | 0.00 | 22,792 |
| 6482 | 168,551 | 44.1 | 45,138 | 2,875 | 7,390 | 796 | 1,054 | 13,593 | **1,183** | **989** | 2.19 | 18,247 |
| 6493 | 161,178 | 57.4 | 53,132 | 4,811 | 6,570 | 1,121 | 902 | 18,640 | **396** | **650** | 1.22 | 20,692 |
| 6521 | 479,006 | 42.6 | 109,538 | 6,353 | 15,635 | 882 | 1,453 | 32,671 | **1,375** | **1,478** | 1.35 | 51,169 |
| 6523 | 600,995 | 52.4 | 157,044 | 7,872 | 25,721 | 444 | 1,586 | 53,643 | **528** | 0 | 0.00 | 67,250 |
| 6533 | 97,237 | 41.7 | 19,838 | 1,004 | 2,261 | 145 | 299 | 5,295 | **88** | **14** | 0.07 | 10,746 |
| 6539 | 218,842 | 51.8 | 57,388 | 3,482 | 15,395 | 1,614 | 2,248 | 14,276 | **452** | 0 | 0.00 | 19,921 |
| 6551 | 174,918 | 39.6 | 40,096 | 6,515 | 8,688 | 801 | 3,616 | 10,198 | **352** | **617** | 1.54 | 9,926 |
| 6566 | 692,878 | 52.6 | 244,050 | 7,708 | 28,846 | 3,075 | 2,121 | 61,145 | **488** | **3,242** | 1.33 | 140,667 |
| 6579 | 577,283 | 49.4 | 147,871 | 8,986 | 18,695 | 2,413 | 2,036 | 41,399 | **956** | 0 | 0.00 | 73,386 |
| 6580 | 119,600 | 30.4 | 14,436 | 251 | 721 | 31 | 64 | 2,178 | **31** | **275** | 1.90 | 11,160 |
| 6593 | 793,534 | 58.8 | 254,564 | 56,084 | 46,752 | 26,725 | 16,534 | 72,951 | **1,852** | **956** | 0.38 | 33,666 |
| **20** | **6,007,789** | | **1,659,255** | **155,421** | **267,589** | **53,762** | **51,012** | **474,798** | **14,609** | **27,019** | | **642,064** |

**`class_NK_ILC = 0` is gone.** The per-cell gate names NK/ILC cells on **all 20 donors**; the
sub-cluster level names an NK/ILC sub-cluster on **13 of 20**. Five donors that v2 reported as
`Immune_NK_ILC = 0` — 6380, 6493, 6523, 6533, 6551 — carry 369, 396, 528, 88 and 352 individually
corroborated NK/ILC cells, and three of them carry a named NK sub-cluster as well.

**Where the NK cells were hiding.** Of the 14,609 per-cell NK/ILC calls, **6,110 sit in clusters
this build's own stage 2 called `Immune_T`** and 8,499 in clusters it called `Immune_NK_ILC`
(`controls.csv`, `nk_calls_from_Immune_T_clusters`). On **6374, 6380, 6523 and 6533 every single NK
call comes from a T-called cluster** — there is no NK-called cluster on those sections at all, which
is exactly the shape of the defect: at the cluster unit those cells had nowhere to go.

**Internal coherence, asserted at write time.** The seven T_NK labels sum to the compartment on
every donor, at both units; the producer raises rather than writes if they do not. The rest of the
release sums too: 3,102,878 Immune + 1,967,990 Non_immune + 811,330 Uncalled_ambiguous + 45,119 not
partitioned + 80,472 no-nucleus = 6,007,789. **A lymph node being mostly immune is an expectation
and is not a check** — 6580 at 30.4% and 6593 at 58.8% are reported as measured.

### The two units disagree, and that is the measurement

`T_NK_unresolved` is 38.7% of the compartment at the per-cell level and 54.2% at the sub-cluster
level. They are not the same residue and they are not competing estimates of one number:

- The per-cell gate **undercounts by construction**. A CD3-negative cell needs two NK transcripts
  present to be called, and most NK cells at this depth have none.
- The sub-cluster call **cannot see an NK cell that is not part of an NK-enriched sub-cluster**. A
  node's NK cells are scattered through the T zone, and a scattered population does not form a
  sub-cluster.
- **4,838 of the 14,609 per-cell NK calls (33.1%) fall inside a sub-cluster the other rule also
  called NK.** Where a sub-cluster is strongly named, the two agree hard: `6442_t13` is 2,747 cells,
  13 NK genes standing up, and 1,389 of its cells clear the per-cell gate on their own.

---

## 3. What the calls rest on: the partition, rebuilt without the no-nucleus cells (D3)

v2 consumed the shared clustering track's partition and then dropped `nucleus_count == 0` cells from
the per-cluster means. Those cells were inside the Leiden subsample and inside the k-NN vote, so they
had already shaped the clusters the calls rest on. **A cell with no nucleus is excluded from every
calculation, and that is not satisfied by excluding it from the last one.**

`--stage partition` re-runs the clustering track's own algorithm — `cluster_pln.leiden_on_subsample`,
imported rather than copied, same seed derivation, same QC — on the nucleated cells only, into this
release's own tree.

| | cells |
|---|---:|
| in the h5ads | 6,007,789 |
| passing the partition QC (`nCount_Xenium >= 10 and nFeature_Xenium >= 5`) — **what v2's partition contained** | 5,962,253 |
| of those, `nucleus_count == 0` — **removed from the partition here** | **80,055** |
| **partitioned, v3** | **5,882,198** |
| labelled `Unavailable_no_nucleus` in the release (includes 417 that also failed QC) | 80,472 |
| `no_nucleus_cells_inside_partition`, asserted per section | **0** |

The brief put this at 41,685. Measured now it is **80,055**, because **all 22 sections have a local
vendor bundle today and v2 could only reach 14** — the other eight were recorded UNMEASURABLE, never
as passed, and the bundles have since landed. 6436 alone contributes 14,758 and 6593 11,930. The
producer refuses to call a section whose partition still contains a no-nucleus cell.

---

## 4. Detection: one predicate, and what it is measured to cost (D2)

```
lognorm(c,g) = log1p(raw(c,g) * 1e4 / depth(c))      detected  <=>  lognorm >= log1p(1)
```

The target is **fixed at 1e4**, not the section median. This is the **only** detection predicate in
the lane: every detection fraction in the release — per cell, per cluster, per sub-cluster, in the
controls — comes from it, and the producer asserts the constant at run time.

The defect was that two different predicates wore the word "detection" in one release: the clustering
track this lane consumed reports its `det__` columns at `lognorm > 0` under a **section-median**
normalisation. v3 derives every detection fraction it reports from its own predicate and **measures
the divergence rather than describing it**:

- **0 of 6,007,789 cells have depth above 10,000**, so the predicate is exactly `raw >= 1` on every
  cell in this cohort.
- **0 of 121,466,945 raw-nonzero entries in the named gene set fall below it.** `> 0` and
  `>= log1p(1)` are measured to coincide, entry for entry, on this data. The divergence is zero and
  now it is a number rather than an assumption.

*(That count is reported with a float64 accumulator. A float32 sum over 6593's 27.4 million stored
entries is not exact above 2²⁴ and returned a count one higher than the entries it was counting.)*

### The Y-probe control on that predicate

`TBL1Y`, `USP9Y`, `EIF1AY` are the floor group: zero-template in a female donor and, unlike `DDX3Y`
and `KDM5D`, not cross-hybridising with a highly expressed X paralog. In a female donor every
detection of them is false.

**Nine female donors, eleven sections** — the brief said two; measured, it is nine.

| section | donor | floor group | all five Y probes |
|---|---|---:|---:|
| 6579b | 6579 | **0.0027%** | 0.0287% |
| HDL115 | 115 | 0.0029% | 0.0279% |
| 6533 | 6533 | 0.0052% | 0.0293% |
| 6476 | 6476 | 0.0063% | 0.0499% |
| 6523 | 6523 | 0.0064% | 0.0295% |
| 6579c | 6579 | 0.0065% | 0.0398% |
| HDL112 | 112 | 0.0069% | 0.0380% |
| 6380 | 6380 | 0.0073% | 0.0270% |
| 6374 | 6374 | 0.0148% | 0.0959% |
| 6482 | 6482 | 0.0173% | 0.0836% |
| 6579a | 6579 | **0.0174%** | 0.0894% |

**0.0027%–0.0174% against the ~0.3% requirement — inside it by a factor of 17 to 111.** In the
eleven male sections the same probes run 3.49%–23.17%, which is template, not background; for male
donors this quantity is **UNMEASURED** and the release says so rather than writing 0.

### The control the negative gate needs, and it is not flattering everywhere

A negative call is only as good as the probe's background. Both ends are measured per section, on
that section's own cells, in `controls.csv`:

- **CD3 detection in a compartment that is definitionally CD3-negative** (B/Plasma) is the gate's
  false-positive floor: **4.0% (6461) to 35.0% (6593)**. B cells do not express CD3; whatever this
  rate is, it is background.
- **CD3 detection in T_NK** is its ceiling: **22.5% (6580) to 86.1% (6593)**.
- Where those two are close the per-cell CD3 negative is weak. On 6593 the floor is 35.0% against a
  ceiling of 86.1%; on 6461 it is 4.0% against 41.5%. **The per-cell NK count on 6593 should be read
  with that in mind**, and it is the reason the sub-cluster estimate exists.
- **The NK rule applied verbatim to a compartment that is not NK** (CD3-negative cells of the
  B/Plasma compartment with ≥2 NK naming genes) runs **0.07%–1.79%**, against **0.28%–17.35%** in
  T_NK. The NK evidence is 3–10× enriched where it should be, on every section.

---

## 5. The NK/ILC sub-clusters, gene by gene

All 14, with the NK-exclusive naming genes that stood up in each and the two independent
corroborations — how often that sub-cluster's own cells detect any NK naming gene against the rest of
the compartment, and how many of its cells clear the per-cell gate.

| sub-cluster | cells | NK genes up | CD3 det. (rest of T_NK) | NK-set det. (rest) | cells clearing the per-cell gate |
|---|---:|---:|---|---|---:|
| `6436_t00` | 10,078 | 14 | 65.1% (67.7%) | 49.3% (12.8%) | 464 |
| `HDL112_t06` | 3,600 | 17 | 26.3% (61.4%) | 80.4% (15.4%) | 1,214 |
| `6566_t22` | 3,242 | 16 | 34.9% (39.6%) | 22.8% (5.5%) | 147 |
| `6442_t13` | 2,747 | 13 | 35.7% (82.7%) | 94.8% (30.0%) | 1,389 |
| `6521_t07` | 1,478 | 5 | 26.7% (48.5%) | 41.8% (19.6%) | 107 |
| `6461_t06` | 1,060 | 17 | 17.1% (38.7%) | 82.7% (14.7%) | 456 |
| `6482_t03` | 989 | 15 | 46.3% (53.5%) | 93.8% (23.7%) | 379 |
| `6593_t04` | 956 | 17 | 69.5% (83.7%) | 95.6% (34.8%) | 220 |
| `6347_t06` | 842 | 16 | 41.8% (62.8%) | 72.6% (20.0%) | 216 |
| `6493_t05` | 650 | 8 | 47.2% (55.8%) | 52.0% (14.9%) | 51 |
| `6551_t10` | 617 | 6 | 61.6% (72.5%) | 68.9% (24.9%) | 56 |
| `6442_t11` | 471 | 5 | 54.6% (81.1%) | 93.6% (32.1%) | 139 |
| `6580_t19` | 275 | 7 | 4.4% (22.1%) | 6.9% (5.7%) | 0 |
| `6533_t08` | 14 | 2 | 42.9% (41.8%) | 14.3% (8.5%) | 0 |

`HDL112_t06` stands up on **17 of the 19** NK-exclusive naming genes — `CD244 EOMES FCGR3A KIR2DL1
KIR2DL4 KIR3DL1 KIT KLRB1 KLRC1 KLRC2 KLRD1 KLRF1 NCAM1 NCR1 NCR3 RORC TBX21` — with CD3 detection at
26.3% against 61.4% in the rest of its compartment. That is not one probe doing the work; it is the
core.

**Two are thin and are flagged, not hidden.** `6533_t08` is 14 cells on `KIT` and `KLRB1` alone and
carries `nk_call_thin_two_or_three_genes`. `6580_t19` has seven genes but only 6.9% of its cells
detect any of them, on a median depth of 34 — that is the lowest-depth section in the cohort and the
call is correspondingly weak.

**Removed by the set-level detection contrast (§1) — 9 sub-clusters, 28,361 cells.** They are named
here because the earlier rule counted them:

| sub-cluster | cells | NK genes up | NK-set det. | rest of T_NK |
|---|---:|---:|---:|---:|
| `6579a_t01` | 9,276 | 2 (`KIR2DL1 KIR3DL1`) | 7.2% | 14.6% |
| `HDL115_t05` | 7,339 | 2 (`KIT NCR3`) | 6.5% | 7.7% |
| `6579b_t08` | 2,897 | 2 (`KIR3DL1 NCAM1`) | 3.1% | 10.3% |
| `6374_t07` | 2,656 | 3 | 13.3% | 18.2% |
| `6566_t23` | 2,499 | 3 | 3.7% | 5.7% |
| `6436_t02` | 1,855 | 3 | 11.6% | 15.6% |
| `6580_t11` | 695 | 3 | 4.0% | 5.8% |
| `6580_t04` | 661 | 2 | 4.5% | 5.7% |
| `6580_t15` | 483 | 2 | 3.9% | 5.7% |

Every one has cells carrying **less** NK evidence than the compartment average. That is the whole
reason the contrast is there, and it is why 115, 6374, 6476, 6523, 6539 and 6579 have no NK
sub-cluster: not "no NK cells" — their per-cell gate finds 179 to 956 each — but **no sub-cluster on
those sections separates them.**

---

## 6. The coarse stage, unchanged in method

Stage 1 (immune vs non-immune on clusters) and stage 2 (identity-core argmax within immune, admitted
only when a gene that may name the winning lineage stands out) are v2's, run on the new partition.
**362 clusters, 22 sections.**

| stage 1 | clusters | cells |
|---|---:|---:|
| Immune | 140 | 3,102,878 |
| Non_immune | 159 | 1,967,990 |
| Uncalled_ambiguous | 63 | 811,330 |

| stage 2 | clusters | cells |
|---|---:|---:|
| Immune_T | 46 | 1,354,750 |
| Immune_B_Plasma | 33 | 978,063 |
| Immune_NK_ILC | 19 | 304,505 |
| Myeloid | 33 | 413,344 |
| Granulocyte | 8 | 34,200 |
| Immune_uncalled | 1 | 18,016 |

`Immune_T` and `Immune_NK_ILC` together are the **T_NK compartment**; their cells never keep those
labels in the per-cell release — they carry the stage-3 split, with the cluster's own call preserved
beside it in `cluster_lineage`.

**Against v2, on the 6,007,859 `cell_id` joins both releases cover** (70 more than the 6,007,789
cells, because 6579's 35 duplicated ids join twice)**:** stage 1 agrees on 4,945,612 (82.3%), the
coarse cluster lineage on 4,650,331 (77.4%). The partition changed, so this is a comparison of two
partitions as much as of two rules; agreement is not evidence for either and is stated so the size of
the change is on the record. Cross-tabs per donor are in `compare_v2/`. **v2's 208,491
`Immune_NK_ILC` cells were a cluster-level statement, and v2's own verification measured CD3E at
59.5% detection across them** — the reason this build splits that compartment rather than counting
it.

---

## 7. Cross-modality, compositional (§ `cross_modality/`)

**18 donors have a pLN section on both platforms; 11 have a typed one on both today.** Xenium-only:
112, 6347, 6374, 6533, 6551, 6566, 6579, 6580, 6593. PhenoCycler-only: 6450, 6505, 6534.

These are **different sections of different blocks**. The join is on `donor_id` and nothing else — no
cell matched to a cell, no coordinate, no spatial claim — and the comparison is on **fractions**,
never counts, because the two assays sample different areas. **Neither side was tuned to move the
agreement and nothing here was accepted because it agreed.**

PhenoCycler source: `data/development/typing/pln/pheno/per_donor_counts.csv`, mtime
`2026-09-17T20:03:16Z`, sha256 `6bc3b62e…`, columns `tnk_T_CD4 tnk_T_CD8 tnk_T_reg tnk_NK_ILC
tnk_T_NK_unresolved`. That lane is moving today; the release records the exact file this ran against.

| quantity | Pearson | Spearman | Xenium median | PhenoCycler median | median \|diff\| |
|---|---:|---:|---:|---:|---:|
| % immune of all cells | 0.552 | 0.464 | 52.4% | 69.0% | 10.9 pts |
| T_NK share of immune | 0.489 | 0.400 | 0.531 | 0.277 | 0.299 |
| NK/ILC share of T_NK — Xenium per-cell | −0.485 | −0.400 | 0.0079 | 0.0167 | 0.026 |
| **NK/ILC share of T_NK — Xenium sub-cluster** | **−0.624** | **−0.722** | **0.0122** | **0.0167** | **0.046** |
| CD8 share of CD4+CD8 | 0.549 | 0.627 | 0.694 | 0.465 | 0.202 |
| Treg share of resolved T | 0.219 | 0.218 | 0.072 | 0.231 | 0.115 |

### The disagreements, which are the interesting part

1. **Xenium calls fewer cells immune, on all 11 donors.** Median 52.4% against 69.0%, up to 32.9
   points apart on 6461 (41.0% vs 73.9%). Eleven of eleven in the same direction, so it is a
   property of the two methods, not noise.
2. **Xenium puts a much larger share of the immune compartment in T_NK** — 0.531 vs 0.277 median, on
   every one of the 11 donors. The two sides are not dividing the immune compartment the same way.
3. **Xenium's CD8 share of CD4+CD8 is higher on all 11 donors** (median 0.694 vs 0.465), but it is
   the **best-correlated** subset quantity in the table (Spearman 0.627). The level differs
   systematically; the donor ordering survives.
4. **Xenium's Treg share is lower on 10 of 11** (median 0.072 vs 0.231) and the correlation is near
   zero (0.219). `FOXP3` is a low-abundance transcript and `FOXP3` protein is a nuclear stain; these
   are not measuring the same thing with the same sensitivity, and the release does not pretend they
   are.
5. **NK/ILC: the levels land close and the ordering is anti-correlated.** Medians 0.0122 (Xenium
   sub-cluster) against 0.0167 (PhenoCycler) — both around 1–2% of T_NK — but the per-donor
   Spearman is **−0.72**, and the pattern behind it is almost exact:

   | | Xenium sub-cluster | PhenoCycler |
   |---|---:|---:|
   | 6436 · 6461 · 6442 · 6482 | 7.5% · 5.4% · 4.6% · 2.2% | **0 on all four** |
   | 115 · 6380 · 6476 · 6523 · 6539 | **0 on all five** | 2.3% · 1.7% · 2.9% · 5.6% · 5.5% |
   | 6493 · 6521 | 1.2% · 1.4% | 1.2% · 6.9% |

   The two sides find an NK population on **almost disjoint sets of donors**. On the nine donors at
   the top of that table one side or the other is at zero, and a zero on either side is a
   sub-cluster-level statement — *this section's partition did not separate an NK population* — not
   a measured absence. Xenium's own per-cell gate finds 179 to 528 individually corroborated NK
   cells on each of its five zeros, so it is not that the cells are absent. **This row is a
   measurement of two partitions, and it should be recomputed when the PhenoCycler lane re-runs.**
   It is reported as it stands, not withheld and not smoothed, and neither side was moved to
   improve it.

The two lanes share no threshold, no gene, no marker and no line of code. That is why the comparison
is worth making, and also why a difference cannot be charged to one side without a further
measurement.

---

## 8. A partial run cannot land on a cohort-wide name (D4)

`build_pln_immune_pheno.py --donors X --stage standout` rewrote three cohort-wide files with the one
donor it had processed. v2 had the same shape: `--only 6533` rewrote `release.json`, the cluster
tables and the summaries for that subset alone, and §9 of the previous version of this document said
so in prose instead of fixing it.

Every cohort-shaped write in v3 goes through `partial_run_guard.RunScope`. A run whose section set is
a strict subset of what the lane could process writes to
`_partial/sections=<slug>/` with a `SCOPE.json` beside it, and the cohort file is not touched, not
read and not merged. **What the lane could process is every section an eligible donor has an h5ad
for — not what happens to be partitioned today**; scoping against the latter would let a one-section
run call itself whole. Verified: `--sections 6347 --stage call` announced the partial scope and wrote
21 files under `_partial/sections=6347/`.

---

## 9. Extension: the donors now transferring

Sections are **discovered** from `/data/xenium/h5ad_out/*_pLNXenium.h5ad` and donor eligibility comes
from `phenocycler.cohort.filter_eligible_donors` (`6457` refused there; this file does not restate
the list). A section already partitioned is skipped unless `--force`, and the cohort rollups are
rebuilt from whatever per-section artifacts are on disk.

Adding a donor is: drop the h5ad in, add its bundle row to `xenium_paths.csv`, re-run the four
stages. Nothing in the producer changes. `6450`, `6505`, `6534`, `6591` and `6623` have no pLN h5ad
on this machine; nothing about them is reported as zero.

---

## 10. What this does not license

- **No spatial claim** of any kind, and none across the pancreas/pLN boundary. Nothing here computes
  a distance, adjacency or neighbourhood, and the cross-modality comparison is compositional.
- **No cell-to-cell matching across modalities.** §7 joins on `donor_id` alone.
- **No non-immune lineage is resolved.** `Non_immune` is terminal; the four non-immune cores appear
  only as far as stage 1 needs them to exclude.
- **Nothing claims acinar or ductal.** Exocrine is one label in this project and would mean nothing
  in a node.
- **NK is not split from ILC.** The identity core is one combined core, so the label is combined;
  the ILC-hint detections are carried per cell as description.
- **The per-cell NK count is a floor, and the sub-cluster count is a sub-cluster-level statement.**
  Neither is "the number of NK cells in the node". `n_nk_exclusive_naming_detected` is written per
  cell so a reader can set their own floor.
- **The regional genes are evidence, not calls.** `CR2 CXCL13 CCL19 CCL21 BCL6 AICDA SELL CD27 MKI67`
  are in `clusters/cluster_gene_stats.parquet`; the regional reading is
  `docs/PLN_REGIONAL_CLUSTERING_2026-09-17.md`'s.

---

## 11. Reproduce

```bash
cd /home/smith6jt/IO60panc2nd/Islet-Explorer-Senior
.venv/bin/python analysis/xenium/scripts/pln_immune_typing_v3.py --stage partition   # ~13 min
.venv/bin/python analysis/xenium/scripts/pln_immune_typing_v3.py --stage call        # ~2 min
.venv/bin/python analysis/xenium/scripts/pln_immune_typing_v3.py --stage tnk         # ~5 min
.venv/bin/python analysis/xenium/scripts/pln_immune_typing_v3.py --stage cross
.venv/bin/python analysis/xenium/scripts/pln_immune_typing_v3.py --stage compare
# one section:  --sections 6533   (writes under _partial/, never on a cohort name)
```

CPU, repo `.venv`, no GPU. `partition` and `tnk` take a random subsample and record the seed they
derived; `call`, `cross` and `compare` take no random draw at all.

## 12. Files

| path | |
|---|---|
| `analysis/xenium/scripts/pln_immune_typing_v3.py` | **the producer** |
| `…/release_pln_immune_v3/release.json` | method, panel, controls, nucleus rule, 22 input sha256 |
| `…/release_pln_immune_v3/label_space.json` | every label with the reason it exists |
| `…/release_pln_immune_v3/calls/donor_id=<d>/roi=pln_1/data.parquet` | per cell |
| `…/release_pln_immune_v3/partition/` | the nucleus-clean partition and its per-section reports |
| `…/release_pln_immune_v3/tnk/` | T_NK sub-cluster assignments and the 249 sub-cluster profiles |
| `…/release_pln_immune_v3/per_donor_summary.csv` · `per_donor_tnk.csv` · `controls.csv` | |
| `…/release_pln_immune_v3/cross_modality/` | §7, with the PhenoCycler file's sha256 and mtime |
| `…/release_pln_immune_v3/compare_v2/` | §6 cross-tabs |
| `data/development/xenium_typing/release_pln_immune_v2/**` | superseded, **left in place** |
| `data/development/typing/pln/clusters/xenium/**` | the shared clustering track — **read only** |
| `core/panels/data/identity_core/*.csv` | the one panel authority; no copy written |

**Keying.** `cell_id` is the Xenium barcode verbatim, `-1` suffix included, within `donor_id`/`roi`.
One measured exception: **6579's three sections share 35 `cell_id` values**, recorded in
`release.json` at `duplicate_cell_ids_within_donor`; for that donor the key is
(`section_key`, `cell_id`). No other donor has any.
