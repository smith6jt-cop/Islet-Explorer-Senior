# Unsupervised regional clustering of the pLN, both modalities — 2026-09-17

Builder: `core/phenocycler/scripts/pln/cluster_pln.py`
Outputs: `data/development/typing/pln/clusters/{phenocycler,xenium}/` (gitignored, 345 MB)

This is a **description of the measurements**, not a typing and not a release. It does not call a
cell positive, does not compute a reference level, and nothing downstream is gated, tuned or
validated by any cluster in it. Existing QuPath classifications, disease labels and prior calls were
never read into a feature matrix.

Every number below was produced by running the builder and is reproduced from a file it wrote; the
file is named beside the number.

---

## 1. "Tissue does not mean anything at all" — measured, and what follows from it

In the pancreas, QuPath's pixel-classifier objects hand you the context for free: islet, blood
vessel, lymphatic, nerve, ordinary tissue. A reference level is then fitted *inside* a context, and
it means "above the ordinary tissue this structure sits in".

**A pLN section has no such partition.** In `core/phenocycler/data/integration/cells_pln`:

| | measured |
|---|---|
| columns | 128, and **none of them is `context`** |
| `parent_raw` | one distinct value, `'Annotation (Tissue)'`, on **all 8,961,168 rows**, 0 nulls |
| `qupath_classification` | **NULL on all 8,961,168 rows** |

"Tissue" here names the whole section, so it names nothing. There is no ordinary tissue to draw a
reference population from because there is no partition of the section at all — not a coarse one,
not a degraded one, none.

**So clustering is the substitute for context.** A lymph node's regional cell types are distinct
enough that an unsupervised partition *discovers* the structure the pancreas is handed. The naming
that follows is then internal to one section: a cluster is named by how its canonical markers stand
against **the other clusters on that same section**. No level, quantile, threshold or scale factor
is imported from anywhere — not from the pancreas, not from another slide, not from a pooled cohort.
Section 4 is that table and Section 4.1 is its contract.

### Withdrawn, and not reintroduced anywhere in this build

- **Any number of pancreas origin reaching a pLN call** — level, threshold, scale factor, quantile.
  The earlier per-donor pancreas-p99 transfer is withdrawn; its numbers are a property of that
  mistake and are not findings, and none of them is used or quoted here as one.
- **`subpanel.fit_levels` / `reference_level` / `context == 'tissue'` applied to the pLN.** An
  earlier version of this document proposed the opposite of the pancreas transfer — fitting the
  levels on the pLN slide's own reference population. That is withdrawn too: it synthesises the
  reference population the table above says is absent. Neither direction survives.
- **`IMMUNE_REQUIRE_MIN = 0.10` as an absolute floor on pLN values.** It is defined against a
  reference level that does not exist here.

## 2. What was clustered

Within donor, always. Never jointly across modalities, never a shared embedding, never a co-plot in
space, and no comparison of pLN to pancreas in space — different blocks on one slide.

| | PhenoCycler pLN | Xenium pLN |
|---|---|---|
| source | `core/phenocycler/data/integration/cells_pln` | `/data/xenium/h5ad_out/*_pLNXenium.h5ad` |
| in the root | 8,961,168 rows, 18 donors | 24 h5ads |
| **clustered** | **5,122,269 cells, 14 donors, 293 clusters** | **5,962,253 cells, 20 donors / 22 sections, 364 clusters** |
| key | `object_id` (QuPath UUID), preserved verbatim | `cell_id`, preserved verbatim; join key is (`section_key`, `cell_id`) |
| features | 54 markers (56 minus DAPI and Beta_actin) | top-2000 HVG **within each section** |
| transform | log1p → within-donor robust z (median / 1.4826·MAD), clipped [−5, 10] | normalize_total → log1p → scale(max 10) |
| clusters per unit | 17–24 | 11–39 |
| k-NN holdout agreement | 0.930 min, 0.939 median, 0.950 max | 0.840 min, 0.920 median, 0.948 max |

Source: `cluster_report.json` and `prose_facts.json` in each modality's directory.

**Eligibility, and why an object count is not a cell count.** `analysis_reason` in `cells_pln`
already encodes it, measured over all 8,961,168 rows: `''` 5,122,269 (analysis-eligible, the ones
clustered) · `no_measurements_exported` 2,786,845 · `missing_nucleus_geometry` 518,429 ·
`cell_area_below_min` 509,390 · `cell_area_above_max` 14,207 · `low_solidity` 10,028. The root was
built from 10,192,986 QuPath objects; **that is an object count and is never reported as a cell
count.** A cell with no nucleus is excluded from every calculation here, which is the 518,429.

**Who is missing, and why.** Four PhenoCycler pLN donors — **6551, 6566, 6593, 6623** — have zero
intensity measurements in the export, so they cannot be clustered. That is an absent measurement,
**not zero immune cells**, and they need re-running in QuPath. Three of the four (6551, 6566, 6593)
are covered by the Xenium lane, so only 6623 is absent from both. 6579 is absent from the
PhenoCycler pLN root because `pln_node_annotations.csv` still records `excluded=True` for it while
`cohort.py` reinstated it on 2026-09-17; it is present in Xenium as three sections (6579a/b/c). 6457
is refused by `phenocycler.cohort` in both.

Donors clustered in **both** modalities: 11 (115, 6380, 6436, 6442, 6461, 6476, 6482, 6493, 6521,
6523, 6539).

**Measurement handling.** Zero PhenoCycler cells were dropped for a partial measurement, and no
marker was unavailable or zero-MAD on any of the 14 slides, so all 54 features were used on every
donor. In Xenium, 45,536 of 6,007,789 cells were dropped for `nCount_Xenium < 10` or fewer than 5
genes detected. Nothing was imputed; a missing measurement would have removed the feature or the
cell, never become a zero. **The Xenium h5ads carry no nucleus geometry**, so the no-nucleus
exclusion has no per-cell handle there — recorded as *unavailable*, not as passed.

## 3. Determinism, run rather than claimed

The previous seeding was `SEED + hash(donor)`. CPython salts `str` hashing per process unless
`PYTHONHASHSEED` is set, so the documented commands did not reproduce their own output. It is now a
BLAKE2b digest of the fixed seed and the donor/section string (`stable_seed`). A second, quieter
instance of the same defect was fixed with it: `eta_squared` summed its levels in `set` order, so a
saturated eta came back as `1.0` in one interpreter and `1.0000000000000002` in the next, which
moved a permutation p. Levels are now summed in sorted order and the result is clipped to [0, 1].

`--stage determinism` runs the same command in two fresh interpreters **with deliberately different
`PYTHONHASHSEED` values** — the exact variable the old seeding was riding on — and compares the
assignment parquet bytes. Measured (`determinism_check.json`):

| modality | units | sha256 identical across the two interpreters | reproduces the tracked build |
|---|---|---|---|
| PhenoCycler | 115 (212,390 cells), 6461 (176,411), 6534 (124,782) | **yes, all three** | **yes, all three** |
| Xenium | 6347 (79,361 cells), 6380 (261,960) | **yes, both** | **yes, both** |

"Reproduces the tracked build" is a join on `object_id` / `cell_id` between the fresh run and the
tracked assignments: every row joined and every `cluster_id` matched.

**This is also why every number in this document was regenerated.** The previous revision quoted
292 PhenoCycler and 362 Xenium clusters, 18–24 and 11–38 per unit, a 0.938 / 0.917 median holdout,
and a marker table topped by Caveolin 8.87 / CD68 8.83. None of those is what the current files say.
They came from a run made before the seeding was fixed — which is precisely the failure the fix
addresses, showing up in the prose.

## 4. What the clusters express, relative to the other clusters on the same section

### 4.1 The table the naming track consumes

`naming_profiles.parquet` (293 and 364 rows) and `naming_profiles/<unit>.parquet`, built by
`--stage profile` from the stored assignments. Per marker/gene per cluster:

| column | what it is |
|---|---|
| `med__` | median in the cluster, on the measured scale (PhenoCycler raw intensity; Xenium log-normalized) |
| `mad__`, `iqr__` | 1.4826·MAD and p75−p25 **within** the cluster — the heterogeneity handle |
| `logmed__` (PhenoCycler) / `mean__`, `det__` (Xenium) | the central value the contrast is taken on, and Xenium detection fraction |
| **`vsrest__`** | **the naming contrast**: this cluster's central value minus the median of the **other clusters' central values on the same section** |
| `sepz__` | `vsrest__` divided by 1.4826·MAD of those other clusters' central values |
| `othermad__` | `sepz__`'s denominator, written beside it |
| `n_cells`, `frac_of_section` | cluster size |

Every quantity is a comparison between clusters measured in the same channel on the same section.
That is what makes it calibration-free across slides, and it is the whole point: nothing has to mean
the same thing on two slides.

**`sepz__` is never read alone, and `othermad__` exists to make that checkable.** A ratio with a
near-zero denominator is large for a reason that is not separation. Measured: **AICDA reaches
`sepz` 49.4 on a `vsrest` of +0.011** log-normalized units, because a gene detected in 0.08% of
cells has essentially no between-cluster spread to divide by. Read the two together.

### 4.2 PhenoCycler — how far a marker separates one cluster from its section's others

From `naming_profiles.parquet`, over the 293 clusters. `max vsrest` is on the log1p intensity scale;
`clusters sepz>3` counts clusters and `donors` counts the donors they come from.

Canonical immune markers (`lineages.LINEAGE_MARKERS` minus `SUPPORTING_ONLY`):

| marker | max vsrest | max sepz | othermad there | clusters sepz>3 | donors |
|---|---|---|---|---|---|
| CD68 | +4.00 | 19.9 | 0.194 | 30 | 13 |
| CD209 | +2.18 | 12.1 | 0.149 | 40 | 14 |
| CD163 | +3.75 | 11.6 | 0.194 | 25 | 11 |
| CD66 | +1.59 | 11.2 | 0.137 | 11 | 10 |
| MPO | +4.10 | 10.2 | 0.388 | 15 | 12 |
| CD206 | +2.26 | 9.9 | 0.201 | 31 | 11 |
| CD20 | +1.93 | 8.3 | 0.231 | 12 | 9 |
| Iba1 | +2.37 | 5.3 | 0.347 | 13 | 8 |
| CD79a | +2.16 | 5.1 | 0.364 | 9 | 6 |
| CD11c | +2.11 | 4.5 | 0.466 | 8 | 5 |
| **CD3e** | **+3.09** | **3.7** | 0.830 | **3** | **1** |
| **CD4** | +1.38 | 3.7 | 0.290 | **3** | **1** |
| **CD8** | +3.84 | **2.4** | 1.619 | **0** | **0** |

Structural, for exclusion: Caveolin +3.67 / sepz 18.3 (55 clusters, 14 donors), CD34 +2.47 / 10.9
(20, 14), B3TUBB +3.74 / 8.8 (3, 3), Collagen_IV +1.88 / 7.9, CD31 +2.99 / 7.4, Podoplanin +2.43 /
5.9, Pan_Cytokeratin +1.60 / 5.3, SMA +4.88 / 3.7, EpCAM and E_cadherin reach no cluster at sepz>3.

Regional and state: Ki67 +3.26 / 15.6 (17 clusters, 13 donors), CD21 +2.51 / 6.9 (9, 7), HLA_DR
+2.30 / 6.3, Granzyme_B +1.71 / 6.3, PCNA +1.70 / 5.8, PD_1 +1.91 / 4.2 (1, 1), TCF_1 +1.56 / 3.8
(3, 2); TOX, CD57, ICOS, Bcl_2 and LAG3 reach no cluster at sepz>3.

**The T markers are the weakest separators on this panel, and the between-cluster contrast says so
more gently than the within-donor z does.** Against the donor's own median cell
(`marker_separation.csv`), no cluster in any of the 14 donors reaches z > 2 on CD3e (max 1.49), CD4
(1.36), CD8 (1.10) or TCF_1 (1.82), and TOX clears it in exactly one cluster of one donor (2.26). Against the section's other clusters, CD3e does reach +3.09 and
sepz 3.7, but in **one** donor and three clusters, and CD8 clears sepz 3 nowhere despite a +3.84
vsrest, because its between-cluster spread (othermad 1.62) is the largest on the panel.

This is what a node does to a T-marker measurement: the ordinary tissue *is* lymphocytes, so the
contrast between a T-rich cluster and the rest of the section is small even when the cells are
there. **That is a statement about how much separation the measurement offers, not about how many T
cells the node contains, and it is not evidence for importing a level from anywhere.** It is the
reason the naming track must read CD3e/CD4/CD8 together with their `othermad__`, and the reason a
T-zone cluster in these data is recognised by a *combination* rather than by one marker clearing a
line.

### 4.3 Xenium — the named canonical and regional genes

From `naming_profiles.parquet` over the 364 clusters; `vsrest` in log-normalized units. `clusters
vsrest>0.25` and `sections` count where a substantial contrast recurs.

| gene | max vsrest | max sepz | othermad there | clusters vsrest>0.25 | sections (of 22) | best cluster detection |
|---|---|---|---|---|---|---|
| CCL21 | +2.10 | 37.9 | 0.0385 | 43 | **21** | 98.5% |
| CR2 | +1.67 | 92.2 | 0.0096 | 8 | 4 | 95.2% |
| MS4A1 | +1.65 | 90.9 | 0.0062 | 33 | **21** | 97.5% |
| NKG7 | +1.52 | 191.1 | 0.0038 | 7 | 6 | 95.3% |
| CD19 | +1.30 | 112.6 | 0.0038 | 28 | 20 | 93.3% |
| CD79A | +1.20 | 105.9 | 0.0044 | 39 | 19 | 93.4% |
| CCL19 | +1.06 | 105.1 | 0.0100 | 12 | 11 | 85.4% |
| CD3E | +1.02 | 12.6 | 0.0348 | 48 | 18 | 88.9% |
| CXCL13 | +0.88 | 217.3 | 0.0041 | 4 | 4 | 87.2% |
| CD68 | +0.84 | 79.6 | 0.0007 | 16 | 11 | 80.2% |
| CD4 | +0.83 | 24.0 | 0.0228 | 23 | 13 | 89.3% |
| PTPRC | +0.74 | 5.2 | 0.0748 | 55 | 19 | 96.1% |
| CD8A | +0.68 | 47.3 | 0.0100 | 17 | 14 | 67.6% |
| ITGAX | +0.68 | 219.3 | 0.0015 | 10 | 8 | 73.3% |
| MKI67 | +0.36 | 93.3 | 0.0018 | 1 | 1 | 54.4% |
| FOXP3 | +0.32 | 157.5 | 0.0009 | 3 | 3 | 37.1% |
| SELL | +0.29 | 11.0 | 0.0042 | 1 | 1 | 43.6% |
| CD27 | +0.28 | 23.7 | 0.0071 | 2 | 2 | 47.1% |
| BCL6 | +0.25 | 53.9 | 0.0014 | 0 | 0 | 26.3% |
| ITGAM | +0.23 | 72.4 | 0.0031 | 0 | 0 | 39.2% |
| KLRD1 | +0.21 | 88.9 | 0.0023 | 0 | 0 | 29.5% |
| CD3G | +0.16 | 17.1 | 0.0042 | 0 | 0 | 31.5% |
| **AICDA** | **+0.011** | 49.4 | 0.0002 | **0** | **0** | **2.9%** |

Cell-weighted detection over all 5,962,253 clustered cells (`gene_detection_and_enrichment.csv`):
PTPRC 52.1%, CCL21 37.9%, CD3E 28.6%, MS4A1 20.1%, CD79A 18.1%, CD4 14.9%, CD19 14.2%, CD8A 11.8%,
CCL19 9.9%, SELL 9.4%, CD27 8.9%, CD68 7.4%, CR2 6.6%, CD3G 5.1%, NKG7 2.7%, CXCL13 2.7%, BCL6 2.6%,
ITGAX 2.0%, FOXP3 1.9%, ITGAM 1.6%, MKI67 1.1%, KLRD1 0.5%, **AICDA 0.08%**.

**CD3D and LYZ are absent from every one of the 24 pLN h5ads** (measured; the 24 share 5,098 genes,
with three swapping between bundles — CPE/HADH/PCSK2 against B2M/IAPP/MT2A). Nothing was substituted
for them. Every other named gene is present in all 24.

**AICDA is not measurable in these data and BCL6 barely is.** Whatever germinal-centre structure
these nodes contain, this clustering did not resolve it by its two canonical genes, and a
germinal-centre claim from this build would be an assumption rather than a measurement.

**Xenium detection here uses lognorm ≥ log1p(1), not a raw-count predicate.** A raw-count rule of
the `raw >= depth/target` form makes the floor depend on a cell's own depth, and it is what lets
gDNA alone anchor a lineage core — a bias that scales with core size, and the immune cores are the
biggest. `det__` in the profile table is "non-zero after normalize_total → log1p", uniformly.

## 5. The recurring groups

Per-donor clusters were related across donors afterwards by correlating their marker profiles
(average linkage on correlation distance, cut at 14 groups for both modalities — the same value for
each, so the parameter is not tuned per result). This describes which per-donor clusters resemble
each other. No cell was ever placed in a shared space with a cell from another donor.

**PhenoCycler — 5 of the 14 groups appear in all 14 donors**, 4,618,136 cells (90.2% of clustered
cells). From `meta_clusters.csv`:

| group | donors | cells | top features (mean within-donor z) |
|---|---|---|---|
| 1 | 14 | 1,574,862 | CD206 +1.49, CD209 +1.26, CD68 +1.15, Iba1 +0.74, CD163 +0.70, CD14 +0.66 |
| 12 | 14 | 1,358,105 | TCF_1 +0.84, CD3e +0.70, FOXP3 +0.69, PD_1 +0.66, CD56 +0.62, CD4 +0.61 |
| 2 | 14 | 843,159 | Caveolin +4.27, CD34 +2.57, Collagen_IV +1.21, CD31 +1.10, SMA +0.50 |
| 9 | 14 | 701,054 | CD20 +1.75, CD79a +1.52, CD21 +1.35, HLA_DR +1.03, Podoplanin +0.56 |
| 13 | 14 | 140,956 | Ki67 +5.28, PCNA +2.13, Granzyme_B +1.71, CD20 +0.66, CD21 +0.56 |

Group 4 (IDO1 +3.16, HLA_DR +2.98, CD11c +1.98) reaches 13 donors and 97,968 cells; group 6 nine
donors; group 5 eight. Four groups rest on 1–3 donors and are named as such in the table — group 3
is 6539 alone (789 cells, Collagen_IV +7.03), group 7 is 6521 alone, group 14 is 6523 alone (69,017
cells), group 10 is three donors.

**Xenium — 6 of the 14 groups appear in all 20 donors**, 4,972,750 cells (83.4%):

| group | donors | cells | top features (cluster mean − section mean, log-normalized) |
|---|---|---|---|
| 10 | 20 | 1,620,527 | PTPRC +0.21, CD3E +0.19, CD8A +0.11, CD247 +0.09, SEMA4D +0.09, ZAP70 +0.08 |
| 3 | 20 | 988,447 | MS4A1 +0.65, CD19 +0.41, CD79A +0.40, CD22 +0.36, TNFRSF13C +0.32, BANK1 +0.21 |
| 2 | 20 | 956,512 | SHANK3 +0.32, CCL14 +0.30, TGFBR2 +0.28, EPAS1 +0.27, PECAM1 +0.27, PLVAP +0.24 |
| 5 | 20 | 714,992 | **CCL21 +0.91, CCL19 +0.23**, PDGFRB +0.15, MYLK +0.10, S1PR3 +0.07 |
| 6 | 20 | 496,869 | CD14 +0.35, CD163 +0.34, SLCO2B1 +0.23, CD4 +0.18, MRC1 +0.18, CD68 +0.16 |
| 4 | 20 | 195,403 | XBP1 +0.63, MZB1 +0.57, POU2AF1 +0.30, PDIA4 +0.25, DERL3 +0.22, CD38 +0.16 |

Group 9 (SOST/CPVL/CLEC9A/ITGAX, 680,589 cells) reaches 19 donors and group 1
(PDGFRB/MYLK/MCAM/NOTCH3, 294,845) reaches 18. Six small groups hold 14,069 cells between them;
group 8 (IL1RL1 +0.77, GATA2 +0.73, KIT +0.68, CTSG +0.58) is 4 donors and 6,205 cells and the other
five are one donor each — four of them 6580.

**What is consistent with what, stated separately.** The groups above are named by expression, and
the anatomical reading is an inference: group 5's CCL21/CCL19 dominance is what a T-zone
fibroblastic reticular network would look like; group 3's MS4A1/CD19/CD79A is what a B follicle
would look like; PhenoCycler group 9's CD20/CD79a/CD21 is the same statement in protein. **None of
those names was assigned to a cluster** — clusters are `<unit>_c<NN>` and groups are numbers — and
none of them was used to accept or reject anything.

## 6. Relation to the node layer

`/data/islet_layers/pln` — measured node outlines over 18 donors, joined by `object_id`. The join
covers all 5,122,269 clustered PhenoCycler cells: **every clustered cell falls inside a node on its
slide** (`node_relation.csv`).

The Xenium clusters have **no relation to this layer and none was invented**. The layer is a
PhenoCycler-frame geometry (`frame = phenocycler_um:<donor>:pln_1`) and does not address Xenium
cells, which come from separately segmented bundles with a disjoint id namespace.

Node counts among the 14 clustered donors: 7 donors have 3 nodes, 4 have 2, and **115 and 6380 have
exactly one**, so a per-node statement is vacuous for them. **6534 is excluded from every per-node
comparison**: its layer records `node_scope = 'slide'` — its one surviving outline is the section —
so a per-node number for 6534 is a per-slide number wearing a node's name.

That leaves **234 donor-clusters across 11 multi-node donors** as testable.

- **1 of the 234 is confined to a single node**: 6442_c19, **53 cells**, in `pLN_top`. One cluster
  of 53 cells is a slide effect if it is anything; it is named rather than rounded to zero.
- median maximum log2(observed / expected-by-node-size) = **+0.57** (1.5×)
- 115 of 234 clusters reach 1.5× in some node; **36 of 234 reach 2×**; the largest is **+3.07**
  (8.4×)
- median total variation distance between a cluster's node distribution and the node-size
  distribution = **0.204**

Read plainly: clusters are shared across the nodes of a donor, with a modest and common tilt toward
one of them, and a thin tail of strongly node-skewed small clusters. Nothing here is a
node-specific population — but the maintainer should note this is also the shape you would see if
two outlines on a slide are sampling the same node twice.

## 7. Batch eta-squared on every cross-donor quantity

**The batch variable is a proxy and that must be read with the numbers.**
`data/donor_acquisition_runs.csv` and `data/donor_acquisition_batch.csv` record the **pancreas**
scans (`6442_Scan1.er.qptiff`). Most pLN slides are separate images (`6442pLN_Scan1.er.qptiff`) that
no tracked manifest in this checkout records. A few donors are the exception — their recorded image
is a combined block — and `covariate_provenance.csv` counts them, measured from the manifests rather
than asserted: **1 of the 14 PhenoCycler donors** (115, `HDL115pancLN_Scan1.er.qptiff`) and **3 of
the 20 Xenium donors** (112, 115, 6533). For those the covariate *is* the pLN acquisition. For the
other 13 and 17 it is not, and 6347 and 6580 have no image record at all. Everywhere else
`stain_group` / `acq_run` / `coarse_batch` / `software` are **donor-level grouping variables, not a
verified pLN scan batch.** A high value is still a warning; a low one is not a clearance.

**Raw eta-squared is unreadable in a cohort this size and is not reported alone.** With 14 donors
over 9 stain groups and 11 acquisition dates, a near-donor-unique grouping partitions essentially
all the variance whatever the values are. Every row of `batch_eta_squared.csv` therefore carries the
median of 2,000 random relabellings with the same level sizes, and a permutation p.

17 cross-donor quantities per modality: `composition_PC1..PC3` and the per-donor fraction of each of
the 14 recurring groups.

| modality | grouping | levels | donors used | max eta² | null median | min p | quantities p<0.05 |
|---|---|---|---|---|---|---|---|
| PhenoCycler | stain_group | 9 | 14 | 1.000 | 0.615 | 0.032 | 2 of 17 |
| PhenoCycler | acq_run | 11 | 14 | 1.000 | 0.802 | 0.023 | 1 of 17 |
| PhenoCycler | coarse_batch | 2 | 14 | 0.210 | 0.038 | 0.095 | **0 of 17** |
| PhenoCycler | software | 3 | 14 | 0.432 | 0.123 | 0.040 | 2 of 17 |
| PhenoCycler | disease_status | 3 | 14 | 0.487 | 0.119 | 0.022 | 1 of 17 |
| Xenium | stain_group | 8 | 16 | 0.792 | 0.461 | 0.022 | 1 of 17 |
| Xenium | acq_run | 11 | 16 | 0.999 | 0.684 | 0.052 | **0 of 17** |
| Xenium | coarse_batch | 3 | 18 | 0.514 | 0.085 | 0.002 | 1 of 17 |
| Xenium | software | 3 | 16 | 0.640 | 0.100 | 0.001 | 1 of 17 |
| Xenium | disease_status | 3 | 18 | 0.819 | 0.087 | 0.012 | 2 of 17 |

Donors with no recorded value for a grouping are **dropped from that grouping**, never pooled into
an "unrecorded" level — that would invent a level made of donors sharing only an absent record.
Xenium loses 4 donors on `stain_group` / `acq_run` / `software` (112, 6347, 6533, 6580 — not in
`donor_acquisition_runs.csv`) and 2 on `coarse_batch` / `disease_status` (6347, 6580).

**The composition axes, stated plainly rather than called unstructured.** The previous build's
`composition_PC2` sat at eta² 0.975 on `acq_run` while its write-up said the composition was not
batch-structured. That was wrong to write, and this build measures the same axis afresh
(`batch_eta_squared.csv`):

| quantity | grouping | eta² | null median | p |
|---|---|---|---|---|
| PhenoCycler composition_PC1 | acq_run | **0.937** | 0.795 | 0.110 |
| PhenoCycler composition_PC2 | acq_run | **0.895** | 0.779 | 0.202 |
| PhenoCycler composition_PC1 | stain_group | **0.838** | 0.630 | 0.082 |
| PhenoCycler composition_PC2 | stain_group | **0.830** | 0.615 | 0.082 |
| PhenoCycler composition_PC2 | software | **0.432** | 0.113 | **0.046** |
| PhenoCycler composition_PC1 | disease_status | 0.368 | 0.119 | 0.077 |
| Xenium composition_PC2 | stain_group | 0.523 | 0.469 | 0.385 |
| Xenium composition_PC2 | acq_run | 0.673 | 0.686 | 0.533 |

**The PhenoCycler composition carries eta² 0.84–0.94 on the acquisition variables.** Those values do
not clear their own null — an 11-level grouping over 14 donors returns a median 0.795 on random
labels — so they are not evidence of a batch axis. **They are also not a clearance**, and with 14
donors this design cannot separate the two. `composition_PC2` on `software` (0.432, p 0.046) is the
one composition axis that clears, and a composition axis tracking the acquisition software is not a
biological axis. The honest position: the composition is reported, it is not corrected, and it
**gates nothing**.

**Quantities that clear p<0.05** (`batch_eta_squared.csv`, `prose_facts.json → batch.significant`):

*PhenoCycler.* group 1 fraction on `acq_run` (eta² 0.980, null 0.812, p 0.023) and on `stain_group`
(0.904, null 0.638, p 0.046) — group 1 is the CD206/CD209/CD68 group, the largest in the build;
group 4 fraction on `stain_group` (0.890, null 0.617, p 0.032); group 11 fraction on `software`
(0.364, p 0.040); `composition_PC2` on `software` (0.432, p 0.046); group 2 fraction on
`disease_status` (0.487, null 0.119, p 0.022). **A myeloid-fraction difference that tracks
acquisition date and stain group should not be read as a myeloid result**, and group 2 is the
Caveolin/CD34 structural group, where a disease signal and a staining difference look identical.

*Xenium.* group 6 (CD14/CD163/MRC1/CD68) carries `software` 0.640 p 0.001, `coarse_batch` 0.514
p 0.002 and `stain_group` 0.792 p 0.022 — three flags on one quantity, and the recorded covariates
describe the *pancreas* acquisition, not the Xenium run. group 8 carries the highest disease eta² in
the build (0.819, p 0.012) but rests on **4 donors and 6,205 cells**; group 1 carries
`disease_status` 0.506 p 0.039. `coarse_batch` and `disease_status` are the alias CLAUDE.md warns
about and cannot be separated here.

**No disease claim is made from any of these.** Nothing was dropped for carrying a confound; every
confounded quantity is reported with it, and none of them gates anything.

### 7.1 The batch term on the quantities the two immune lanes PUBLISH

Section 7 above measures the batch term on the *clustering's own* cross-donor quantities. It was
absent from the quantity a reader actually meets: **the per-donor immune numbers each lane
publishes.** Those were reported per donor with no batch term beside them. They now have one, for
every quantity in both published tables, computed by
`cluster_pln.py --modality <m> --stage published-eta` into
`published_quantity_batch_eta_squared.csv` / `.json` under each modality's cluster root.

The inputs are the tables the lanes publish, read and not rewritten:
`data/development/typing/pln/pheno/per_donor_counts.csv` (14 donors, 16 quantities) and
`data/development/xenium_typing/release_pln_immune_v2/per_donor_summary.csv` (20 donors,
21 quantities). Every quantity is a **fraction** of a stated denominator, because eta-squared on a
count would partly measure how many cells each node happened to yield. Both tables were checked for
internal coherence before anything was computed on them, and both hold: PhenoCycler
immune + non-immune + uncalled = eligible, the five `class_*` columns sum to immune, and the
published `pct_immune` equals its own ratio; Xenium immune + non-immune + uncalled_ambiguous +
not_clustered + no_nucleus_excluded = cells, and the ten `s2__` columns sum to cells.

**The headline quantity, both lanes.** PhenoCycler `immune_fraction_of_eligible`, Xenium
`immune_fraction_of_cells`:

| modality | grouping | levels | donors | eta² | null median | p |
|---|---|---|---|---|---|---|
| PhenoCycler | acq_run | 11 | 14 | **0.848** | 0.779 | 0.320 |
| PhenoCycler | stain_group | 9 | 14 | **0.764** | 0.610 | 0.174 |
| PhenoCycler | software | 3 | 14 | 0.110 | 0.110 | 0.501 |
| PhenoCycler | coarse_batch | 2 | 14 | 0.084 | 0.034 | 0.300 |
| PhenoCycler | disease_status | 3 | 14 | 0.067 | 0.115 | 0.680 |
| Xenium | acq_run | 11 | 16 | 0.650 | 0.678 | 0.566 |
| Xenium | stain_group | 8 | 16 | 0.499 | 0.459 | 0.406 |
| Xenium | disease_status | 3 | 18 | 0.132 | 0.106 | 0.399 |
| Xenium | software | 3 | 16 | 0.054 | 0.108 | 0.686 |
| Xenium | **xenium_batch** | 2 | 20 | 0.020 | 0.024 | 0.534 |
| Xenium | coarse_batch | 3 | 18 | 0.017 | 0.095 | 0.882 |

**0.848 is elevated and it is not separable from its own null at this n, and that is the finding.**
An 11-level grouping over 14 donors returns a median eta² of 0.779 on labels carrying no structure
at all, so 0.848 is a number those degrees of freedom can produce on their own (p 0.320). The same
holds for `stain_group` at 0.764 against a null of 0.610 (p 0.174). Neither is evidence of an
acquisition effect on the immune fraction and **neither is a clearance**: with 14 donors spread over
11 acquisition dates this design cannot tell the two apart. Reporting 0.848 alone would manufacture
a batch effect out of degrees of freedom; calling it "not significant" and moving on would hide one.

**The Xenium lane had no published batch term at all until now**, and it is the one lane that can
carry a batch variable measured on the node itself: `xenium_batch` is the instrument batch of the
pLN section, read from the `roi` rows beginning `pln` in
`core/phenocycler/data/integration/xenium_paths.csv` (2 levels — 6KTRAF and Y9QRJ6 — over all 20
donors, with no donor dropped). It is the only grouping in this document that is not a pancreas-scan
proxy. The immune fraction returns eta² 0.020 against a null of 0.024 (p 0.534).

**What clears p<0.05, stated rather than filtered out.** Nothing below is dropped, corrected, or used
to gate anything.

| modality | quantity | grouping | donors | levels | eta² | null | p |
|---|---|---|---|---|---|---|---|
| PhenoCycler | non_immune_fraction_of_eligible | stain_group | 14 | 9 | 0.874 | 0.621 | 0.041 |
| PhenoCycler | class_T_cell_fraction_of_immune | stain_group | 14 | 9 | 0.871 | 0.630 | 0.049 |
| PhenoCycler | class_Granulocyte_fraction_of_eligible | coarse_batch | 14 | 2 | 0.318 | 0.043 | 0.045 |
| PhenoCycler | class_Granulocyte_fraction_of_immune | coarse_batch | 14 | 2 | 0.307 | 0.045 | 0.046 |
| Xenium | not_clustered_fraction_of_cells | software | 16 | 3 | 0.698 | 0.095 | 0.0005 |
| Xenium | not_clustered_fraction_of_cells | coarse_batch | 18 | 3 | 0.652 | 0.049 | 0.0045 |
| Xenium | not_clustered_fraction_of_cells | stain_group | 16 | 8 | 0.925 | 0.382 | 0.049 |
| Xenium | s2__Immune_NK_ILC_fraction_of_cells | stain_group | 16 | 8 | 0.772 | 0.436 | 0.041 |
| Xenium | s2__Immune_NK_ILC_fraction_of_immune | stain_group | 16 | 8 | 0.750 | 0.458 | 0.038 |
| Xenium | s2__Immune_B_Plasma_fraction_of_immune | **xenium_batch** | 20 | 2 | 0.262 | 0.025 | 0.023 |

Four things follow, none of them a verdict on a marker or a lane:

- **The granulocyte fraction is the PhenoCycler quantity most tied to the coarse batch**, on a
  2-level grouping whose null is 0.043. `coarse_batch` is aliased with disease status at exactly
  this level, so that row cannot be read as either. It sits directly beside defect **D1** — CD66 is
  on the builder's structural-veto list while `lineages.py` names it a granulocyte canonical, and
  6523 returns 0 granulocytes against 1,334–68,985 elsewhere. A quantity with a known defect in its
  producer and a batch flag is not a result.
- **`not_clustered_fraction_of_cells` carries the strongest flags in the whole table** (p 0.0005 on
  `software`, 0.0045 on `coarse_batch`). It is not a biological quantity — it counts cells the
  partition did not reach — and it tracking the acquisition is what one would expect of a
  technical quantity. It is published because a denominator that moves with batch moves every
  fraction computed on it.
- **The Xenium NK/ILC fraction flags on `stain_group`** (0.772 / 0.750, p ≈ 0.04, against a null of
  0.44–0.46). The NK/ILC quantity is the one under active revision by the nk-t-subsets track; this
  row is a reason to read its per-donor spread with the batch term visible, not a reason to
  withhold it.
- **`s2__Immune_B_Plasma_fraction_of_immune` is the only quantity to clear on `xenium_batch`**
  (0.262, null 0.025, p 0.023) — the only grouping here measured on the pLN section itself. Two
  levels, 20 donors, nothing dropped. It is the single cleanest batch signal in this document and it
  belongs to the B/plasma fraction, not to the immune fraction.

The PhenoCycler `class_NK_ILC` column is 0 on all 14 donors in the table read here, so every
`class_NK_ILC` row comes back with zero variance and its eta-squared is `nan` — **unavailable
because the quantity is constant, not because the batch term is absent.** The rows are in the CSV
with `state` recording it. When the nk-t-subsets track lands the CD3e−/CD8+ split, rerunning
`--stage published-eta` fills them in with no change to this lane.

## 8. Coherence checks

Internal consistency only. None of these passes by reproducing a biological finding, and none would
fail because a lymph node contained the wrong thing (`naming_profiles_coherence.json`):

| check | PhenoCycler | Xenium |
|---|---|---|
| assigned rows | 5,122,269 | 5,962,253 |
| key unique over the modality | **true** (`object_id`) | false — `cell_id` repeats across bundles |
| key unique within the unit | true | **true**; join key is (`section_key`, `cell_id`) |
| cluster ids in assignments = cluster ids in the profile table | true | true |
| max disagreement between a cluster's size in the assignments and in the profile table | **0** | **0** |
| profile sizes sum | 5,122,269 | 5,962,253 |
| max abs(Σ frac_of_section − 1) over units | 0.0 | 0.0 |
| assigned cells that are not analysis-eligible in `cells_pln` | **0** of 5,122,269 | n/a |
| max abs(median over a section's clusters of `vsrest__`) | 0.372 (6523, SMA) | 0.076 (6551, PTPRC) |

The last row is an indicator, not a law: `vsrest` is a value minus the median of the others, so its
median over a section's clusters lands near zero unless that section's clusters are strongly
asymmetric on that marker. 6523's SMA is that case, and the measured value is a property of the
data, not a tolerance that was met.

## 9. What this build does not license

- **No spatial claim across tissues.** Nothing here compares pLN to pancreas in space, and nothing
  co-plots the two modalities. Pancreas and pLN are different blocks on one slide; the only
  pancreas/pLN comparison the method allows is compositional, joined on `donor_id` and disease
  status.
- **No cluster is an anatomy.** The names are `<unit>_c<NN>` and `meta_cluster <n>`. Section 5
  states what is consistent with what, separately and as inference.
- **No germinal-centre result.** AICDA is at 0.08% detection with a max `vsrest` of +0.011, and BCL6
  clears +0.25 in no cluster of the 364.
- **No disease result.** See section 7.
- **No exocrine or acinar or ductal anything.** Those labels have no meaning in a node, and the
  60-plex could not carry them anyway.
- **This is not the typing.** Naming the clusters by canonical marker/gene enrichment against the
  section's other clusters, then resolving within immune with the `core/panels` identity-core
  subpanels (`03_immune_T_cell`, `04_immune_B_plasma`, `05_immune_NK_ILC`, `06_immune_myeloid`,
  `07_immune_granulocyte`), is the next track. `naming_profiles.parquet` is built for it.
- **The existing `data/development/typing/pln/pheno` calls are not input to anything here.**
  `cluster_immune_link.csv` records where those calls land per cluster, as description. That track's
  own manifest states its stage 1 fitted levels on "the pLN slide's OWN reference population", which
  is withdrawn (section 1), so its agreement with any cluster is not evidence for it and is not read
  as one. It needs rebuilding on the internal contrast before its numbers mean anything.

## 10. Reproducing it

```bash
cd /home/smith6jt/IO60panc2nd/Islet-Explorer-Senior
.venv/bin/python core/phenocycler/scripts/pln/cluster_pln.py --modality phenocycler
.venv/bin/python core/phenocycler/scripts/pln/cluster_pln.py --modality xenium
```

Stages run independently: `--stage cluster`, `profile`, `summarize`, `link`, `prose`,
`published-eta` (section 7.1), and
`--stage determinism`, which is the check in section 3 and writes `determinism_check.json`. Seeded
at 20260917 throughout, derived per unit with `stable_seed` (BLAKE2b), so the commands reproduce
their own output.

### 10.1 A partial run never writes a cohort-wide name

**This bit twice on 2026-09-17.** `build_pln_immune_pheno.py --donors 6539 --stage standout`
rewrote three cohort-wide files with the one donor it had just processed — `cluster_calls.csv`
293 rows → 19, `cluster_heterogeneity.csv` 293 → 19, `per_donor_separation.csv` 14 → 1 — and
earlier the same day `protein_rna --donors 6539` dropped 6436 and 6442 from its summary index the
same way. Nothing errored. The file simply became a smaller file wearing the same name, and the
next reader had no way to tell. `cluster_pln.py --stage profile --donors 115` had the identical
shape: it concatenates the per-donor naming profiles it just wrote and lands them on the
cohort-wide `naming_profiles.parquet`.

**The rule, in `core/phenocycler/scripts/pln/partial_run_guard.py`:** a run whose unit set is a
strict subset of what the lane could process is PARTIAL, and a partial run writes its
cohort-shaped outputs under `_partial/<kind>s=<slug>/`. The cohort-wide path is left exactly as it
was — not read, not merged, not removed.

Three candidates were on the table and the choice is deliberate.

- **Merge** — re-read every per-unit shard on disk and concatenate — was **rejected.** The shards on
  disk are not all from the same generation; a shard from yesterday's parameterisation would merge
  as silently as today's. This repository has already been bitten by a typing tree that stacks
  generations side by side with nothing marking which is which. Merging would turn a *visible*
  truncation into an *invisible* mixture, which is worse.
- **Scope** is what `cluster_pln.py` does. A one-donor naming profile is a useful object; it is just
  not `naming_profiles.parquet`. It is written in full, one directory down, with a `SCOPE.json`
  naming what it covers and what it does not.
- **Refuse** is what `build_pln_layer.py` does, because there is no such thing as a one-donor node
  membership **layer** — a scoped copy of it would be a thing no reader wants. A donor the published
  root already holds and the run does not write stops the build when the reason is a **missing
  input** (its derivation has not been run), and is still removed, and named in the MANIFEST, when
  the reason is a **decision** (`phenocycler.cohort`, or the roster's `excluded` flag). Build a
  partial run into its own root with `--out` instead.

Both are the same rule — the cohort path belongs to the cohort run — and each lane takes the
consequence its artifact deserves.

What is scoped in `cluster_pln.py` is everything that aggregates across donors:
`naming_profiles.parquet`, `naming_profiles_coherence.json`, `cluster_report.json`, `summary.json`,
`regional_enrichment.*`, `cluster_profiles.parquet`, `meta_clusters.csv`,
`meta_cluster_regional_enrichment.csv`, `batch_eta_squared.csv`,
`batch_eta_squared_required_grid.csv`, `donor_meta_cluster_composition.csv`,
`covariate_provenance.csv`, `node_relation.*`, `meta_cluster_node_recurrence.csv`,
`cluster_immune_link.*`, `marker_separation.csv`, `marker_raw_dynamic_range.csv`,
`gene_detection_and_enrichment.csv`, `prose_facts.json`, `determinism_check.json` and
`published_quantity_batch_eta_squared.*`. The per-unit outputs — `assignments/`,
`profiles/<donor>.parquet`, `naming_profiles/<donor>.parquet` — are keyed by their own unit and are
written where they always were.

`summarize`, `link` and `prose` never enumerate their units; they read whatever is on disk. A named
selection is therefore treated as partial for those stages too, because a run that cannot show it
covers the cohort may not write the cohort's file. `--out-root` does **not** make a run whole: a
private root is still a root, and a partial file copied off one must not arrive wearing a cohort
name.

Proven, not asserted. `core/phenocycler/tests/current/test_pln_layer.py` runs a two-donor build and
then a one-donor rebuild into the same populated root and asserts the cohort-wide
`naming_profiles.parquet` is **byte-for-byte** what it was, that the scoped copy holds the one
donor, and that the layer builder refuses when a donor's derivation has gone missing. Both tests
fail against the code as it stood this morning.

Read `cells_pln` with an explicit hive schema — bare `partitioning='hive'` infers `donor_id` as
int32 there, because every pLN donor id happens to be all digits. The builder does this; the tables
it writes are duplicated as parquet for the same reason.

## 11. Open, for the maintainer

1. **6551, 6566, 6593, 6623 need re-running in QuPath.** Zero intensity measurements in the current
   pLN export — 4 of 18 donors and 2,786,845 rows. Absent measurement, not absent immune cells.
2. **6579 disagreement is unresolved.** `cohort.py` reinstated it 2026-09-17; the tracked
   `pln_node_annotations.csv` still records `excluded=True`. The PhenoCycler pLN root withholds it;
   Xenium includes its three sections.
3. **No pLN acquisition manifest exists.** Until one does, every batch statement in section 7 is
   about a donor-level proxy taken from the pancreas scans.
4. **PhenoCycler group 1 — the largest recurring group — tracks acquisition date and stain group**
   (eta² 0.980 p 0.023 and 0.904 p 0.046). It is myeloid-marker-dominated, and that is the pattern a
   staining or exposure difference would also produce.
5. **`class_T_cell` is the wrong name for what the panel resolves, and the rename is not mine
   to make.** The PhenoCycler lane's column is called `class_T_cell`; `lineages.py:31` carries
   the class as `T_NK: (CD3e, CD8, CD4, FOXP3, CD56)`, and `lineages.py:69` already splits it
   `{"CD8_T": (CD8,), "CD4_T": (CD4,), "NK": (CD56,)}`, and what the combination actually resolves
   is CD3e+ CD4+ → CD4 T, CD3e+ CD8+ → CD8 T, CD3e+ FOXP3+ → Treg, **CD3e− CD8+ / CD3e− CD56+ →
   NK / ILC**. A column named `class_T_cell` that carries NK cells in it will be read as a T-cell
   count by everything downstream. The name should be `class_T_NK` at the coarse level, with the
   subsets the nk-t-subsets track is adding named under it — `class_CD4_T`, `class_CD8_T`,
   `class_Treg`, `class_NK_ILC` — so that the subsets sum to the coarse class and the coarse class
   is not a lie about its own contents. **That track owns
   `core/phenocycler/scripts/typing/build_pln_immune_pheno.py`; this is a naming request to it,
   not an edit made here**, to avoid two agents writing the same builder in the same hour. The one
   consumer in this lane is `--stage published-eta`, which reads whatever `class_*` columns the
   table has and names its rows after them, so it follows the rename with no change. The same note
   applies to `class_NK_ILC` being 0 on all 14 donors today: section 7.1 records it as a constant
   quantity rather than as a measured zero batch term.
6. **The pLN typing track's stated stage 1 is withdrawn.** `data/development/typing/pln/pheno` was
   built by fitting levels on the pLN slide's own reference population; section 1 measures why there
   is no such population. That track, not this one, has to change.
