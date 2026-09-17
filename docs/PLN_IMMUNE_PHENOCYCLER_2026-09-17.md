# PhenoCycler pLN — immune populations, named on the clusters

**2026-09-17.** Development lane. Nothing here is a release and the `NO_VERIFIED_TYPING` latch has
not moved.

Builder `core/phenocycler/scripts/typing/build_pln_immune_pheno.py` ·
tests `core/phenocycler/tests/current/test_pln_immune_pheno.py` (58, all passing) ·
output `data/development/typing/pln/pheno/`.

**Three stages.** Stage 1 admits a cluster to the immune compartment, stage 2 names it from the
`core/panels` subpanels, and **stage 3 splits the combined `T_NK` compartment into T_CD4, T_CD8,
T_reg and NK_ILC by negative gating on CD3e**. NK/ILC is named: **25,368 cells on 8 of 14 donors.**

```bash
.venv/bin/python core/phenocycler/scripts/typing/build_pln_immune_pheno.py
PYTHONDONTWRITEBYTECODE=1 .venv/bin/python -m pytest -p no:cacheprovider \
  core/phenocycler/tests/current/test_pln_immune_pheno.py -q
```

## Why this lane has no reference level at all

`cells_pln` has **no `context` column**. `parent_raw` is the single constant string
`Annotation (Tissue)` on all 8,961,168 rows and `qupath_classification` is NULL on every row. The
pancreas's context machinery — islet, blood vessel, lymphatic, nerve, "tissue", all of it read off
QuPath pixel-classifier objects — has nothing to attach to in a lymph node. There is no ordinary
tissue to draw a reference population from because **there is no partition of the section at all**.

So the clustering is the substitute for context. In the pancreas a QuPath object hands you the
context for free and a level is fitted against it; in a node the clustering discovers the structure,
because the regional cell types are distinct. The pipeline is:

1. cluster within donor, per section — `core/phenocycler/scripts/pln/cluster_pln.py`, **consumed
   here, not re-run**;
2. name each cluster by canonical enrichment **relative to the other clusters on that same
   section**;
3. resolve within immune with the five `core/panels/data/identity_core` subpanels;
4. cells inherit their cluster's call, and per-cluster heterogeneity is reported.

## Withdrawn, and gone

| withdrawn | what it was | why |
|---|---|---|
| (a) pancreas → pLN level transfer | each donor's pancreas p99 carried onto their pLN slide | the transferred CD66 level left 1.9 % of the slide beneath it and on 12 donors sat below that slide's own 5th percentile |
| (b) `subpanel.fit_levels` / `reference_level` / `context == 'tissue'` | the level then fitted on the pLN itself, as p99 of "context == tissue" | on this export that population is the whole slide, so **every marker called 1 % of the slide positive by construction** |
| (c) `IMMUNE_REQUIRE_MIN = 0.10` as an absolute floor | positivity in the pancreas lane's scaled space | it is defined against a reference level that does not exist here |

`pancreas_levels()`, `acquisition_factors()`, `transfer_check()` and their
`_reference/*.parquet` outputs were already absent from this checkout before this work; nothing
under `data/` was deleted to keep them so, and the test asserts they stay absent.
`fit_pln_levels()`, `in_fires_units()`, `pln_context()` and `stage1_decide()` — generation (b) —
are now gone from the module, and `subpanel_io` is **not imported at all**.

**Generation (b)'s outputs are still on disk** — `_reference/pln_reference_levels.{parquet,csv}`
and `_reference/reference_population.csv` — because the hard rule for this work was to delete
nothing under any `data/` directory. They are named as superseded in
`_reference/SUPERSEDED_2026-09-17.json`, which also lists the current outputs. **Do not cite them.**
Removing them needs a maintainer decision; leaving a superseded generation beside a current one is
exactly the hazard `typing-tree-stacks-generations` records.

The withdrawn attempts' numbers are a property of the mistake, not of the data, and are not cited
as findings anywhere.

## The statistic: what "stands out" means

For cluster *c* and marker *m*, over the analysis-eligible cells of **one section and that section
only**:

> **AUC[c, m]** = P( a random cell in *c* reads higher on *m* than a random cell elsewhere on this
> same section )

the Mann–Whitney common-language effect size, from that section's own ranks. It was chosen because
it is a comparison *between clusters on one section* (so nothing is imported and nothing has to
mean the same thing on two slides), it is **threshold-free and unit-free** (it reads ranks, so no
intensity value, level, quantile or scale factor enters, and 0.5 is the null), it is robust to the
heavy tails an intensity distribution carries, and its null is *calculable* rather than assumed —
under random relabelling of that section's cells the rank sum of a cluster of size *n* has a known
mean and a tie-corrected variance.

A test asserts the consequence directly: multiplying a channel by 1000 and adding an offset changes
not one call.

**ENRICHED** = AUC > 0.5 **and** tie-corrected rank-sum *p* < 0.05 Bonferroni-corrected over every
cluster × marker test on that section **and** S > 0, where

> **S[c, m] = AUC[c, m] − median over all 56 channels of AUC[c, ·]**

The third conjunct is there because the offset is measured. On donor 115 the per-cluster median AUC
across the panel runs **0.107 to 0.794** — a cluster-wide offset of 0.687 AUC — and a cluster's
EpCAM AUC correlates **0.973** with that median (CD3e, 0.810). With 200k cells a rank-sum *p* is
significant almost for free, so S is what makes "enriched" mean anything. The offset is estimated
over the whole panel, so it cannot be steered by which markers this lane cares about; per-donor
offset spread runs 0.570 (6380) to 0.867 (6505) and is in `per_donor_separation.csv`.

**What S does not do**, stated because it would be easy to claim otherwise: it does not change the
immune-vs-structural arbitration or the subpanel argmax. Both are comparisons between two markers of
the *same* cluster, and a per-cluster constant cancels in a comparison — the calls and the
separation margins are identical with and without it, measured. S changes which markers are admitted
as evidence, nothing else.

## Stage 1 — immune vs non-immune, on clusters

The marker sets are projections of `lineages.LINEAGE_MARKERS` minus `lineages.SUPPORTING_ONLY`, not
hand-picked. **CD45/PTPRC is not on this panel**, so stage 1 cannot be one gate; it is a union of
lineage canonicals arbitrated against the structural set.

| role | markers |
|---|---|
| T_NK | CD3e, CD8, CD4 |
| B_Plasma | CD20, CD79a |
| Myeloid | CD68, Iba1, CD163, CD206, CD11c, CD209 |
| Granulocyte | MPO, CD66 |
| supporting, never alone | CD11b, CD56, FOXP3, M2Gal3 |
| structural, for exclusion | B3TUBB, CD31, CD34, Caveolin, Collagen_IV, E_cadherin, EpCAM, Pan_Cytokeratin, Podoplanin, SMA |

Dropped from structural with the reason travelling with the drop: **Vimentin** (every lymphocyte
carries it; `lineages` puts it in Stromal, which is a statement about the pancreas — excluding on it
would exclude the node) and **TP63** (already `SUPPORTING_ONLY`, so it may not exclude a cluster
alone). `Ker8_18` and `Keratin_5` are wanted by the taxonomy and are not on this export —
unavailable, never zero.

A cluster is **immune** when a canonical immune marker is enriched and reads higher, on S, than any
enriched structural marker; **non-immune** when the structural evidence wins; **uncalled** when
nothing stands out. Uncalled is a real category. In the pancreas such a cell falls through to
Exocrine because exocrine is the background there; a node has no such default, and inventing one
would manufacture a majority class out of nothing. This lane names **no** non-immune lineage — the
structural markers arbitrate and never become a label.

Across 293 clusters: 178 immune, 106 non-immune, 9 uncalled at stage 1.

## Stage 2 — within immune, the subpanels

Membership is `core/panels/data/identity_core/*.csv` reached through the verified crosswalk in
`phenocycler.integration.vocab`, and nothing is patched in by hand — which costs this lane two
canonicals and the loss is recorded rather than papered over. **CD209** is a myeloid canonical in
`lineages` and `06_immune_myeloid.csv` does not list it, so it carries no stage-2 evidence.
**CD66**'s genes are CEACAM1/5/6 and `07_immune_granulocyte.csv` lists CEACAM3, so the granulocyte
subpanel's only naming marker here is **MPO**.

| subpanel | naming markers on this panel |
|---|---|
| T_cell | CD3e, CD4, CD8 |
| B_Plasma | CD20, CD79a |
| NK_ILC | *(none — see below)* |
| Myeloid | CD11c, CD163, CD206, CD68, Iba1 |
| Granulocyte | MPO |

The winner is the strongest subpanel **whatever its enrichment**, and is then required to carry an
enriched canonical of its own lineage. A winner without one is **UNCALLED with that named reason**
and is never reassigned to the runner-up. That is the failure the Xenium draft made — 55.4 % of its
Granulocyte calls detected none of the granulocyte canonicals — and the requirement is written so
it cannot recur. It bound once here: one cluster on 6461, 10,155 cells, `winning subpanel Myeloid
has no enriched canonical of its own lineage`.

**NK_ILC is called on zero clusters AT THIS STAGE, and that zero is a property of the route, not a
measurement.** CD56 is the NK subpanel's only member here and `lineages.SUPPORTING_ONLY` holds it
back from naming a cell alone because it is shared with the Neural lineage. It is therefore not a
stage-1 primary, so a cluster standing out on CD56 alone never enters the immune compartment. The
licence from `subpanel.immune_rules()['Imm-NK-cell']` (CD56 with CD57 or Granzyme B, CD3e negative)
is implemented and tested, and it can only fire for a whole cluster admitted on another lineage's
canonical. No cluster met it.

**That is not where NK is named. Stage 3 names it, from the absence of CD3e — see below.**

A withdrawn generation of this document argued from the zero above that NK could not be named at
all, because no single channel on this panel is NK-specific. **That reasoning is withdrawn.** No
channel has to be. NK cells are CD3e-negative and CD8- and/or CD56-positive; T cells are
CD3e-positive; the panel tells them apart by **combination**, in which absence carries as much as
presence. `lineages` carries `T_NK: (CD3e, CD8, CD4, FOXP3, CD56)` as one class because that is the
**coarse** label, not because the split is out of reach.

## Per donor

All 14 typeable donors separate: every section has clusters carrying an enriched immune canonical,
and no donor had to be withheld for failing to separate.

| donor | eligible | clusters | immune | non-immune | uncalled | % imm | T | B | Myeloid | Granulocyte | Uncalled | max immune S | clusters w/ enriched canonical |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 115 | 212,390 | 19 | 112,876 | 97,293 | 2,221 | 53.1 | 12,224 | 34,899 | 52,691 | 13,062 | 0 | 0.71 | 16/19 |
| 6380 | 219,116 | 19 | 173,106 | 46,010 | 0 | 79.0 | 67,119 | 40,138 | 48,366 | 17,483 | 0 | 0.73 | 16/19 |
| 6436 | 483,690 | 17 | 341,053 | 137,459 | 5,178 | 70.5 | 121,910 | 58,533 | 128,104 | 32,506 | 0 | 0.78 | 14/17 |
| 6442 | 319,242 | 24 | 247,767 | 63,787 | 7,688 | 77.6 | 73,254 | 70,631 | 91,719 | 12,163 | 0 | 0.79 | 22/24 |
| 6450 | 391,651 | 24 | 260,310 | 131,341 | 0 | 66.5 | 79,873 | 54,125 | 93,384 | 32,928 | 0 | 0.80 | 21/24 |
| 6461 | 176,411 | 23 | 130,280 | 46,131 | 0 | 73.9 | 36,056 | 30,358 | 35,101 | 18,610 | 10,155 | 0.83 | 19/23 |
| 6476 | 505,148 | 18 | 344,301 | 160,847 | 0 | 68.2 | 61,586 | 113,885 | 134,015 | 34,815 | 0 | 0.77 | 15/18 |
| 6482 | 325,393 | 23 | 176,390 | 139,792 | 9,211 | 54.2 | 48,782 | 36,946 | 85,973 | 4,689 | 0 | 0.65 | 18/23 |
| 6493 | 221,120 | 20 | 163,236 | 56,735 | 1,149 | 73.8 | 44,659 | 37,463 | 51,931 | 29,183 | 0 | 0.72 | 16/20 |
| 6505 | 236,076 | 21 | 137,360 | 91,431 | 7,285 | 58.2 | 44,022 | 36,512 | 54,657 | 2,169 | 0 | 0.85 | 19/21 |
| 6521 | 748,094 | 24 | 457,937 | 290,157 | 0 | 61.2 | 84,809 | 102,029 | 202,114 | 68,985 | 0 | 0.68 | 19/24 |
| 6523 | 784,013 | 19 | 541,178 | 242,835 | 0 | 69.0 | 189,053 | 121,913 | 230,212 | 0 | 0 | 0.83 | 16/19 |
| 6534 | 124,782 | 21 | 67,431 | 49,825 | 7,526 | 54.0 | 13,994 | 18,137 | 33,966 | 1,334 | 0 | 0.75 | 19/21 |
| 6539 | 375,143 | 21 | 232,995 | 115,551 | 26,597 | 62.1 | 45,685 | 41,387 | 97,545 | 48,378 | 0 | 0.81 | 18/21 |

**Totals: 5,122,269 analysed cells, all clustered — 3,386,220 immune, 1,669,194 non-immune, 66,855
uncalled.** T_cell 923,026 · B_Plasma 796,956 · Myeloid 1,339,778 · Granulocyte 316,305 ·
Uncalled (stage 2) 10,155. The 923,026-cell `T_cell` column is the **combined T_NK compartment**;
stage 3 splits it.

**Never read an object count as a cell count.** The root holds 8,961,168 rows from 10,192,986 QuPath
objects; 5,122,269 are analysis-eligible. `no_measurements_exported` withholds 2,786,845 and
`missing_nucleus_geometry` 518,429 — a cell with no nucleus is excluded from every calculation.

Not typed, each with a measured reason on disk: **6551, 6566, 6593, 6623** — the QuPath export
carries no intensity measurement (7 of 918 fields filled), so there is nothing to cluster and
nothing to name. **Never 0 immune cells.** **6579** is withheld by `pln_node_annotations.csv`.

## Per-cluster heterogeneity

Every cluster carries, in `cluster_heterogeneity.csv`, how concentrated its canonical evidence is at
the cell level. `cell_frac_agrees_class` is the share of a cluster's cells whose own strongest
subpanel is the one the cluster was named; `cell_argmax_entropy` is the normalised Shannon entropy
of that per-cell argmax.

The per-cell comparison uses **one marker per subpanel** — the marker that spoke for that subpanel
on that cluster — because a max over a subpanel's whole set is a max of *k*, and *k* differs: five
naming markers for myeloid against one for granulocyte. Correcting that moved the median granulocyte
agreement from 0.660 to 0.773 and the median T-cell agreement from 0.611 to 0.473, so the
uncorrected statistic was reporting the shape of the panel rather than the mixedness of a cluster.

| class | clusters | cells | median agreement | min | median entropy |
|---|---|---|---|---|---|
| Myeloid | 80 | 1,339,778 | 0.843 | 0.172 | 0.368 |
| Granulocyte | 21 | 316,305 | 0.773 | 0.185 | 0.396 |
| B_Plasma | 44 | 796,956 | 0.716 | 0.308 | 0.552 |
| T_cell | 32 | 923,026 | **0.473** | 0.333 | 0.773 |

**38 of the 177 called clusters — 979,636 cells — have cells that mostly disagree with their own
label**, 18 of them T_cell clusters holding 641,662 cells. The T compartment is the least
concentrated one on this panel: its median cluster agrees with itself on under half its cells. That
is reported, not corrected; a per-cell T call on this lane should be read as "this cell sits in a
cluster the section's CD3e/CD4/CD8 evidence named", not as a per-cell T-cell call.

## Coherence checks — internal contradictions only

No check here compares a result against what a node should contain. A lymph node is expected to be
mostly immune and **that expectation is not a check**; the 53–79 % immune range above is reported,
never verified against.

`_reference/coherence.csv` defines positivity by a **section-internal rank grid** — a cell is
positive for a marker at *q* when it is in the top *q* of that section for that marker — across
q ∈ {0.25, 0.10, 0.05, 0.01}. **No call reads q.** The calls are made on cluster enrichment, so this
is a diagnostic of the diagnostic.

- **The union never passes 100 % on any donor at any q.** CD4-positive *or* CD8-positive exceeding
  the whole T compartment would be an arithmetic impossibility; it does not occur.
- **The sum passes 100 % on 5 of 14 donors, and only at the loosest definition** (q = 0.25):
  6450 124.9 %, 6461 128.7 %, 6482 117.8 %, 6476 112.0 %, 6442 100.2 %, driven by double-positive
  shares of 15.5–40.7 %. At q ≤ 0.10 no donor passes 100 %. The withdrawn attempt's contradiction
  (CD4 69 % / CD8 58 %) is the same shape, and the honest reading is that a definition loose enough
  to call a quarter of a section positive for everything produces it on its own.
- T-cluster cells that are also CD20-positive at q = 0.10: 0.13–5.57 %. B-cluster cells CD3e-positive:
  0.08–8.10 %.
- Immune-cluster cells whose own structural rank beats their own immune rank: 19.3–33.2 %. Stage 1
  arbitrates on the cluster, so this is reported rather than asserted to be zero.

## Where the exclusions came from

AUC is scale-free — it reads ranks and discards magnitude — so a channel that is nearly flat across
a section can still separate a cluster at AUC 0.9, because the ordering of background is still an
ordering. `structural_veto_accounting.csv` therefore reports every exclusion a structural marker
caused **with that marker's measured span on that section** (brightest cluster's median over the
section's median). It gates nothing; turning it into a gate would be an absolute threshold and this
method has none.

| deciding marker | donors | clusters excluded | cells excluded | median span |
|---|---|---|---|---|
| Caveolin | 14 | 30 | 522,378 | 20.3× |
| SMA | 13 | 19 | 280,363 | 5.3× |
| CD34 | 10 | 15 | 229,209 | 8.4× |
| **EpCAM** | 5 | 5 | **222,395** | **1.60×** |
| B3TUBB | 9 | 10 | 155,721 | 2.2× |
| Pan_Cytokeratin | 9 | 13 | 146,227 | 3.2× |
| E_cadherin | 4 | 7 | 60,656 | 2.9× |
| CD31 | 4 | 4 | 42,373 | 10.2× |
| Collagen_IV | 1 | 1 | 7,154 | 1.9× |
| Podoplanin | 2 | 2 | 2,718 | 3.7× |

**Read the EpCAM row before using the non-immune counts.** It decides 222,395 exclusions at a median
span of 1.60× — the narrowest of any deciding channel, an order of magnitude under Caveolin — and on
donor 115 the excluded cluster is the section's largest at 35,251 cells with CD3e at AUC 0.906. The
method cannot indict a channel, only be reported honestly: a rank statistic will arbitrate on a
channel with almost no span exactly as confidently as on one with twenty-fold span, and here it did.
Whether those clusters should have been excluded is a question for the images, not for another
statistic.

## Stage 3 — the T_NK compartment split by negative gating on CD3e

**The deliverable.** The 923,026 cells stage 2 puts in `T_cell` are the combined `T_NK` lineage.
Stage 3 splits them on the presence/absence pattern the panel actually carries:

| pattern | call |
|---|---|
| CD3e+ CD4+ | `T_CD4` |
| CD3e+ CD8+ | `T_CD8` |
| CD3e+ FOXP3+ | `T_reg` |
| **CD3e− CD8+ or CD3e− CD56+** | **`NK_ILC`** |

### The CD3e-negative rule, and why it carries no threshold

Every term is decided on **within-section ranks with an offset removed**, at whichever level the
comparison is made. Nothing imports an intensity, a level, a quantile or a scale factor, and
nothing crosses a section. A test multiplies a channel by 1000 and asserts that not one call moves.

- **Cluster level.** `S = AUC − the cluster's median AUC over the whole panel`, which stages 1 and 2
  already use. **CD3e-negative is `S ≤ 0`:** the cluster reads CD3e no higher than it reads the
  panel at large. `0` is not a cutoff — it is where `S` sits for a marker carrying none of the
  cluster's own level. The raw AUC against the rest of the section is reported beside it.
  *Measured: **not one** of the 32 compartment clusters is CD3e-negative as a whole, by either
  reading.* Every one is a genuine T-cell cluster, so the whole-cluster route never fires and the
  resolution has to be sub-cluster.
- **Cell level.** `r = the cell's place on this section for a marker, minus that cell's own median
  place over the whole panel`. A cell is CD3e-negative when `r_CD3e ≤ 0`, and carries the positive
  NK term when `r_CD8 > 0` or `r_CD56 > 0`.

**The per-cell offset is what makes the gate possible at all, and it is measured, not assumed.**
Inside a cluster a dim cell reads low on CD3e, CD8 and CD56 alike, so a raw comparison compares
brightness rather than identity: on the **raw** ranks the CD3e-low side read *lower* on CD8 in **all
32** compartment clusters (AUC 0.193 to 0.513, never above 0.52), and the gate could never fire.
With each cell's own level removed, CD56 separates on eight of them. It is the same correction stage
1 already applies one level up, and for the same reason.

**A 2-group split is computed and reported, but it decides nothing.** The best 2-group partition of a
cluster maximises between-group sum of squares, which on a largely unimodal cluster cuts it near the
middle: on 6523_c00 it put 44,466 of 101,138 cells on the low side, with *both* sides reading CD3e
above their own panel level (+0.060 and +0.234). NK cells in a T zone are a minority and cannot be
found that way. `split_eta2_cd3e` and `eta2_at_the_zero_cut` are both in the table so a reader can
see whether the zero cut falls anywhere near a natural valley.

**The group-level check is not circular, and it is corrected over its own family.** The two sides are
formed by **CD3e alone**; the check is then whether the CD3e-negative side reads higher than the
CD3e-positive side **of the same cluster** on CD8 or CD56. Only then are that side's CD8/CD56-positive
cells named. The alpha is Bonferroni over the tests *this gate* runs on that section — one per
compartment cluster per NK term, 2 to 8 of them, so α = 0.00625 to 0.025. Borrowing stage 1's
section-wide α (3.7e-5 to 5.3e-5, over 952–1,344 tests) would correct for a family two to three
orders of magnitude larger than the one being run. The count is fixed before any test runs, so no
result can change it.

### Per cluster — what happened, on all 32

`CD8 AUC` and `CD56 AUC` are the CD3e-negative side over the CD3e-positive side of that same cluster,
on offset-corrected ranks. `corr` is the within-cluster Spearman correlation of the two
offset-corrected ranks — reported, never applied.

| donor | cl | cells | CD3e S | CD3e− | of cl | CD8 AUC | CD56 AUC | p | α | corr CD3e/CD56 | named | NK |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 115 | c05 | 12,224 | +0.262 | 302 | 2.5 % | 0.384 | 0.592 | 5e-08 | 0.025 | +0.005 | **yes** | 282 |
| 6380 | c08 | 18,642 | +0.107 | 4,903 | 26.3 % | 0.244 | 0.412 | 1e-75 | 0.00833 | +0.171 | no | 0 |
| 6380 | c11 | 24,462 | +0.169 | 3,949 | 16.1 % | 0.370 | 0.486 | 4e-03 | 0.00833 | +0.061 | no | 0 |
| 6380 | c17 | 24,015 | +0.278 | 1,215 | 5.1 % | 0.314 | 0.531 | 3e-04 | 0.00833 | +0.028 | **yes** | 1,119 |
| 6436 | c11 | 59,752 | +0.159 | 8,114 | 13.6 % | 0.341 | 0.483 | 1e-06 | 0.0125 | +0.054 | no | 0 |
| 6436 | c15 | 62,158 | +0.234 | 3,089 | 5.0 % | 0.244 | 0.379 | 2e-114 | 0.0125 | +0.168 | no | 0 |
| 6442 | c16 | 41,396 | +0.159 | 4,848 | 11.7 % | 0.343 | 0.469 | 2e-12 | 0.0125 | +0.126 | no | 0 |
| 6442 | c17 | 31,858 | +0.291 | 1,102 | 3.5 % | 0.292 | 0.432 | 1e-14 | 0.0125 | +0.105 | no | 0 |
| 6450 | c00 | 43,832 | +0.184 | 5,107 | 11.7 % | 0.240 | 0.427 | 3e-64 | 0.00833 | +0.233 | no | 0 |
| 6450 | c02 | 34,108 | +0.121 | 4,093 | 12.0 % | 0.302 | 0.402 | 5e-92 | 0.00833 | +0.172 | no | 0 |
| 6450 | c11 | 1,933 | +0.085 | 541 | 28.0 % | 0.190 | 0.483 | 3e-01 | 0.00833 | +0.088 | no | 0 |
| 6461 | c11 | 2,023 | +0.321 | 53 | 2.6 % | 0.107 | 0.491 | 8e-01 | 0.00625 | +0.082 | no | 0 |
| 6461 | c15 | 9,641 | +0.197 | 750 | 7.8 % | 0.287 | 0.391 | 3e-23 | 0.00625 | -0.056 | no | 0 |
| 6461 | c16 | 14,875 | +0.207 | 1,194 | 8.0 % | 0.240 | 0.439 | 2e-12 | 0.00625 | +0.116 | no | 0 |
| 6461 | c17 | 9,517 | +0.138 | 258 | 2.7 % | 0.395 | 0.485 | 4e-01 | 0.00625 | +0.279 | no | 0 |
| 6476 | c13 | 61,586 | +0.251 | 1,896 | 3.1 % | 0.343 | 0.516 | 2e-02 | 0.025 | +0.089 | **yes** | 1,808 |
| 6482 | c15 | 9,221 | +0.100 | 2,668 | 28.9 % | 0.198 | 0.472 | 3e-05 | 0.0125 | +0.045 | no | 0 |
| 6482 | c18 | 39,561 | +0.124 | 6,754 | 17.1 % | 0.428 | 0.322 | 1e-78 | 0.0125 | +0.279 | no | 0 |
| 6493 | c10 | 15,437 | +0.294 | 538 | 3.5 % | 0.299 | 0.640 | 3e-28 | 0.0125 | -0.202 | **yes** | 514 |
| 6493 | c13 | 29,222 | +0.238 | 2,192 | 7.5 % | 0.410 | 0.474 | 6e-05 | 0.0125 | -0.016 | no | 0 |
| 6505 | c06 | 16,408 | +0.252 | 770 | 4.7 % | 0.216 | 0.404 | 2e-19 | 0.0125 | +0.162 | no | 0 |
| 6505 | c09 | 27,614 | +0.222 | 825 | 3.0 % | 0.369 | 0.425 | 2e-13 | 0.0125 | +0.190 | no | 0 |
| 6521 | c03 | 35,146 | +0.274 | 2,102 | 6.0 % | 0.202 | 0.435 | 2e-23 | 0.0125 | +0.039 | no | 0 |
| 6521 | c08 | 49,663 | +0.160 | 6,337 | 12.8 % | 0.217 | 0.520 | 3e-07 | 0.0125 | -0.028 | **yes** | 5,816 |
| 6523 | c00 | 101,138 | +0.166 | 11,409 | 11.3 % | 0.357 | 0.566 | 1e-115 | 0.00625 | -0.097 | **yes** | 10,527 |
| 6523 | c01 | 69,017 | +0.108 | 13,271 | 19.2 % | 0.320 | 0.491 | 2e-03 | 0.00625 | +0.038 | no | 0 |
| 6523 | c03 | 15,653 | +0.335 | 116 | 0.7 % | 0.389 | 0.535 | 2e-01 | 0.00625 | -0.055 | no | 0 |
| 6523 | c17 | 3,245 | +0.103 | 740 | 22.8 % | 0.295 | 0.445 | 5e-06 | 0.00625 | +0.053 | no | 0 |
| 6534 | c02 | 11,065 | +0.095 | 2,879 | 26.0 % | 0.239 | 0.518 | 5e-03 | 0.0125 | -0.081 | **yes** | 2,794 |
| 6534 | c13 | 2,929 | +0.129 | 692 | 23.6 % | 0.229 | 0.433 | 1e-07 | 0.0125 | +0.051 | no | 0 |
| 6539 | c06 | 12,079 | +0.147 | 1,970 | 16.3 % | 0.288 | 0.357 | 2e-89 | 0.0125 | +0.200 | no | 0 |
| 6539 | c14 | 33,606 | +0.181 | 2,648 | 7.9 % | 0.332 | 0.590 | 5e-54 | 0.0125 | -0.122 | **yes** | 2,508 |

**CD8 never carries the gate on any of the 32** — `nk_side_auc_cd8` runs 0.107 to 0.428, below 0.5
everywhere, so inside a T cluster CD8 is inseparable from CD3e at the cell level, which is what you
would expect where CD8 T cells dominate the CD8 signal. **CD56 is the term that works**, and where it
does, CD3e and CD56 are disentangled within the cluster: `corr` is +0.089 or below on all eight
clusters that resolve, four of them negative, against a median of +0.097 on the 24 that do not.

**24 clusters refuse, and the refusal carries its numbers.** On those, 59,202 cells do read CD3e at
or below their own panel level *and* CD8 or CD56 above it — and are **not** named NK, because the
side they sit on does not separate from the CD3e-positive side. Each cluster's `resolution_reason`
states the AUC, the p and the α it was held to.

### Per donor — the counts, which sum to the compartment

| donor | T_NK | T_CD4 | T_CD8 | T_reg | NK_ILC | unresolved | CD4 % | CD8 % | CD4+CD8 % |
|---|---|---|---|---|---|---|---|---|---|
| 115 | 12,224 | 6,409 | 4,089 | 1,424 | **282** | 20 | 53.8 % | 34.3 % | 88.1 % |
| 6380 | 67,119 | 24,074 | 26,591 | 15,239 | **1,119** | 96 | 36.5 % | 40.3 % | 76.9 % |
| 6436 | 121,910 | 46,276 | 44,782 | 30,852 | **0** | 0 | 38.0 % | 36.7 % | 74.7 % |
| 6442 | 73,254 | 31,074 | 20,644 | 21,536 | **0** | 0 | 42.4 % | 28.2 % | 70.6 % |
| 6450 | 79,873 | 26,277 | 27,966 | 25,630 | **0** | 0 | 32.9 % | 35.0 % | 67.9 % |
| 6461 | 36,056 | 15,808 | 12,282 | 7,966 | **0** | 0 | 43.8 % | 34.1 % | 77.9 % |
| 6476 | 61,586 | 23,455 | 20,417 | 15,818 | **1,808** | 88 | 39.3 % | 34.2 % | 73.5 % |
| 6482 | 48,782 | 16,317 | 18,480 | 13,985 | **0** | 0 | 33.4 % | 37.9 % | 71.3 % |
| 6493 | 44,659 | 23,294 | 8,727 | 12,100 | **514** | 24 | 52.8 % | 19.8 % | 72.6 % |
| 6505 | 44,022 | 19,686 | 17,056 | 7,280 | **0** | 0 | 44.7 % | 38.7 % | 83.5 % |
| 6521 | 84,809 | 29,709 | 35,657 | 13,106 | **5,816** | 521 | 37.9 % | 45.4 % | 83.3 % |
| 6523 | 189,053 | 81,839 | 75,233 | 20,572 | **10,527** | 882 | 46.1 % | 42.4 % | 88.4 % |
| 6534 | 13,994 | 2,205 | 6,668 | 2,242 | **2,794** | 85 | 19.8 % | 60.0 % | 79.8 % |
| 6539 | 45,685 | 19,173 | 15,602 | 8,262 | **2,508** | 140 | 44.6 % | 36.3 % | 80.8 % |

**Cohort: 923,026 cells — T_CD4 365,596 · T_CD8 334,194 · T_reg 196,012 · NK_ILC 25,368 ·
T_NK_unresolved 1,856.** They sum to 923,026 exactly, on every donor and on the cohort.
**NK/ILC is named on 8 of 14 donors, in 8 clusters.** `T_NK_unresolved` is a CD3e-negative cell in a
resolved cluster that reads neither CD8 nor CD56 above its own level, or a cell carrying no CD3e
measurement at all — an unavailable measurement is never a negative one.

### The coherence check

- **The subsets sum to the compartment.** Every cell of the compartment leaves stage 3 with exactly
  one label, so this is arithmetic; it is asserted per donor in `per_donor_t_nk.csv` and by the test
  suite against the cells on disk.
- **CD4 % + CD8 % within CD3e-positive cells never passes 100 %** — 67.9 % to 88.4 %, because the
  three CD3e-positive terms are an argmax and cannot claim a cell twice. The remainder is `T_reg`.
- Stage 3 labels **no cell outside its compartment**: a myeloid or B cell carries an empty
  `t_nk_class`, checked on disk.

### Six donors return NK_ILC = 0, and that 0 is a measurement

6436, 6442, 6450, 6461, 6482 and 6505. On each the gate **ran** — CD3e, CD8 and CD56 all measured, a
CD3e-negative and a CD3e-positive side both present — and the evidence pointed the other way: their
CD3e-negative cells read *lower* on CD56 (AUC 0.322 to 0.491) as well as on CD8. `nk_ilc_state`
records `measured` for them and `t_nk_resolution.csv` carries the numbers cluster by cluster.

**NK_ILC is written EMPTY, never 0, where it was not decidable** — where the compartment is empty,
where no cluster has both sides with CD8 or CD56 measured, or where every compartment cluster came
back `did_not_resolve`. `nk_ilc_state` then carries the reason and the sentence *"This is NOT zero NK
cells."* No donor is in that state on this cohort; the guard is asserted by
`test_nk_ilc_is_never_written_as_zero_when_the_compartment_did_not_resolve`.

### What stage 3 does not claim

- **`T_reg` is a measurement of FOXP3 on this panel, not a lineage apart from CD4.** A Treg *is* a
  CD4 T cell carrying FOXP3, and `lineages` records FOXP3 above its positivity boundary in 90.4 % of
  cells, which is why it is `SUPPORTING_ONLY` and may not admit a cell. It admits none here: every
  cell reaching the argmax is already inside a cluster admitted on CD3e/CD4/CD8, and `T_reg` is
  written for the CD3e-positive cells that read **highest** on FOXP3. At 21.2 % of the compartment it
  should be read as *the FOXP3-leaning share of the T cells*, not as a Treg count.
  `median_subset_argmax_margin` per cluster (0.050 to 0.409) says how decided that argmax was.
- **`NK_ILC` is one label.** Nothing here separates NK from ILC; the panel carries no channel that
  could. The Xenium identity core `05_immune_NK_ILC.csv` holds 31 genes for it — FCGR3A, KLRD1,
  NKG7, EOMES, CD244, KIR2DL1/4 — and that is where the finer resolution has to come from.
- **No frequency, uniformity or disease ordering entered anywhere.** That eight donors resolve NK and
  six do not is an output, not a check.

## What this lane does not license

- **No spatial claim across tissues.** Pancreas and pLN are different blocks on one slide. Nothing
  here joins them at all — not by distance, adjacency, neighbourhood, KDE or shared scatter, and not
  even by donor. A pancreas/pLN comparison is compositional, joined on `donor_id` and disease status
  and nothing else.
- **No non-immune lineage.** The structural markers exclude; they never become a label.
- **Never acinar or ductal.** The 60-plex carries no marker that could split exocrine, and exocrine
  has no meaning in a node anyway.
- **No result here is a verification target.**

## Outputs

```
data/development/typing/pln/pheno/
  donor_id=<d>/roi=pln_1/data.parquet   per cell: object_id (QuPath UUID, verbatim), cluster_id,
                                        is_immune, immune_class, immune_subtype, call_source,
                                        stage1_rule, stage2_reason, the cluster's standouts, the
                                        cell's own immune/structural section ranks, and stage 3's
                                        t_nk_class / t_nk_resolution_mode /
                                        t_nk_subset_argmax_margin
  donor_id=<d>/roi=pln_1/UNAVAILABLE.json   for a donor that is not typed, with its reason
  cluster_calls.csv                     293 clusters: both stages, every subpanel's standout,
                                        the NK licence, and the raw AUC beside each S
  t_nk_resolution.csv                   the 32 compartment clusters: the CD3e S and AUC, the
                                        CD3e-negative count, both NK terms' AUC and p, the gate's
                                        own alpha, the within-cluster CD3e/CD8 and CD3e/CD56
                                        correlations, the 2-group split that decided nothing, the
                                        mode, the per-cluster subset counts and a prose reason
  per_donor_t_nk.csv                    per donor: the compartment, its four subsets and the
                                        residual, the CD4/CD8 percentages of the CD3e-positive
                                        cells, and nk_ilc_state
  cluster_heterogeneity.csv             per-cluster cell-level agreement and entropy
  per_donor_counts.csv / per_donor_status.csv / per_donor_separation.csv
  structural_veto_accounting.csv
  _reference/coherence.csv, cluster_calls.parquet, SUPERSEDED_2026-09-17.json
  manifest.json
  _partial/<scope>/                     a run that did not cover every clustered donor and every
                                        stage; its PARTIAL_RUN.json names what it covers.  These
                                        are NOT the cohort and nothing above is written from one.
```

## Two defects fixed, 2026-09-17

**D1 — CD66 was listed as a non-immune structural marker.** `lineages.py` defines
`"Granulocyte": ("MPO", "CD66")` and its own docstring records an earlier lane dropping CD66 and
thereby *"folding granulocytes into the myeloid label"*. CD66 nevertheless sat in
`cluster_pln.PHENO_NON_IMMUNE`. It has been removed, and a test now checks **both** that list and
this builder's `STRUCTURAL` set against `lineages` so no immune canonical can return to either.

*Measured effect: **no call in this lane changes**, and that is the finding.* CD66 is **not** in
this builder's structural veto set — `STRUCTURAL` is derived from the four non-immune lineages of
`lineages.py`, and Granulocyte is immune — so no cluster was ever vetoed by it here, and CD66 is a
stage-1 primary that admits 12 clusters, 4 of them as Granulocyte. The stale list was in
`cluster_pln.py`, where it feeds only the **reported** `set__non_immune_structural` column of
`regional_enrichment.csv`, a column nothing in the repository reads. Recomputed from the profiles on
disk, removing CD66 moves that column by up to **0.129** (6450_c16, a Granulocyte cluster, +0.129;
6436_c07 +0.117; 6461_c22 +0.106) and flips its sign on **3 of 293 clusters, 89,442 cells** (6476_c09
17,133 · 6482_c03 28,350 · 6523_c02 43,959). The clustering track has not been re-run — re-running it
would replace the clusters these calls are joined to — so `regional_enrichment.csv` on disk still
carries the old column until it is.

**6523's 0 granulocytes, re-measured: CD66 is not the cause.** Its two MPO-led clusters are excluded
at stage 1 by structural markers — **c05, 92,494 cells**, MPO S 0.065 against EpCAM S 0.082 (a margin
of 0.016), and **c07, 29,866 cells**, MPO S 0.095 against CD34 S 0.536. Both are recorded in
`structural_veto_accounting.csv` with the deciding marker's span on that section. The narrow EpCAM
margin on c05 is worth a look; the CD34 one on c07 is not close. Neither arbitration was changed
here — that would be a different decision from the one this task asked for.

**D4 — a partial run overwrote the cohort-wide tables.** `--donors X` and `--stage standout` wrote
straight over `cluster_calls.csv` (293 rows → 19) and `per_donor_separation.csv` (14 → 1). A run
that does not cover **every** clustered donor **and** every stage now writes to
`_partial/<scope>/` and announces it, leaving the cohort files untouched; it also no longer writes
`UNAVAILABLE.json` markers, which are a cohort-wide statement. Verified by checksum: two partial
runs (`--donors 6442 --stage standout`, then `--stage standout`) left all six cohort tables
byte-identical.
