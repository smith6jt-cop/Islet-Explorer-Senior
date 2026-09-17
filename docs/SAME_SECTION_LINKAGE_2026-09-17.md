# The same-section linkage: how far a cell-level cross-modality statement reaches

**2026-09-17.** Everything numbered here was produced by running
`scripts/senior/same_section_linkage.py` against this checkout. Outputs are in
`data/development/integration/same_section/` (`donor_linkage.csv`, `fanout.csv`,
`agreement.csv`, `linkage.json`).

Everywhere else in this project protein and RNA are **different measurements of different
sections**, and only composition compares. Two lanes claim otherwise. This document establishes
exactly how far each one reaches, and says plainly where it stops.

---

## 1. Two lanes, and they are not the same claim

**The post-Xenium IF re-stain.** After the Xenium run the section is stained for INS and GCG.
The IF image and the transcriptome are therefore the **same physical section**.
`scripts/restain/protein_rna.py` measures the two hormones inside **every Xenium cell's own
polygon** — so on the value axis there is *no segmentation transfer at all*. The values land on
the Xenium cell table as `feat__if_INS` / `feat__if_GCG`. They are **features**: continuous
protein evidence, no positivity asserted, no identity implied. Nobody has defined a true
positive for them (`STAMP_ENDOCRINE_SUBTYPE` is still OFF), and nothing below thresholds them.

**The same-slide lane (6533, 112).** Full PhenoCycler then Xenium on one section, with the
protein *call* carried onto Xenium cells by largest shared area
(`core/phenocycler/scripts/sameslide/lift_calls.py`, read here, not run). That transfer is
many-to-one by construction. A lifted call is **not a measurement of the Xenium cell**; it is a
protein statement about overlapping area, and `lift_calls.py` is explicit about it — it carries
`call_share` and `call_purity` precisely so the reader can see how much the label is worth.

The two lanes fail differently and must not be pooled under one word.

---

## 2. Per donor: what the linkage supports

Gate: a **held-out median nucleus offset** from the refinement's two spatial folds. Above
**2.0 µm** the transform does not support a cell-level statement — a nucleus is roughly 7 µm
across, so at that scale the nearest warped nucleus is routinely a *different* nucleus. This is
a geometric statement about the transform. It is not fitted to any result and no biological
outcome was consulted in setting it.

| donor | supports | held-out p50 / p90 (µm) | IF channels named | IF values stamped |
|---|---|---|---|---|
| **6436** | **per-cell** | 0.614 / 1.623 | yes | `feat__if_INS`, `feat__if_GCG` |
| **6442** | **per-cell** | 0.657 / 1.702 | yes | `feat__if_INS`, `feat__if_GCG` |
| 6539 | per-region | 0.524 / 1.311 | yes | **no — values on disk, never stamped** |
| 6482 | per-region | 0.506 / 1.368 | no | no |
| 6551 | per-region | 0.519 / 1.474 | no | no |
| 6523 | per-region | 0.560 / 1.444 | no | no |
| 6591 | per-region | 0.605 / 1.612 | yes | no |
| 6521 | per-region | 0.672 / 1.794 | no | no |
| 6593 | per-region | 0.746 / 1.827 | no | no |
| 6374 | per-region | 0.875 / 1.805 | no | no |
| 6566 | per-region | 1.859 / 2.877 | no | no |
| 115 | per-region | **2.670** / 4.052 | yes | no |
| 6380, 6461, 6476, 6493 | per-region | **no gate was run** | no | no |
| **112** | **nothing** | — | — | — |
| **6533** | **nothing** | — | — | — |

Counts: **2 per-cell, 14 per-region, 2 nothing.** Of the 18 donors with registration artifacts,
16 have a refined transform and only 10 of those have a held-out gate at all.

### The donors that support nothing at cell level, and why

- **6533** — registration failure. There is no `refined.json`. What survives in
  `data/development/restain/valis/donor_id=6533/` is **4,738 IF nuclei against 1,097,489 Xenium
  nuclei** (1,051,000 cells) — the IF nucleus detection recovered **0.43%** of what the
  transcriptome segmented, so there is almost nothing to fit a transform to. Matched through the
  superseded valis-stage warp, the median nearest-neighbour separation is **311.5 µm**; that
  number is only worth quoting as a diagnosis, since the refined warp it should be measured
  against does not exist. This is not a weak registration, it is the absence of one, and no
  number from this donor may be carried to a cell. Its same-slide acquisition is real and recorded in `section_modes.csv` (PhenoCycler
  70-plex scan `0080870` then Xenium `0080870`), but the lifted calls were produced on
  HiPerGator and are not in this checkout, so nothing here can reproduce or check them.
- **112** — no refined transform and no refined gate, and the same upstream failure: **7,973 IF
  nuclei against 551,335 Xenium nuclei**, or **1.45%**. `if_nuclei_um.npy` and
  `xenium_nuclei_um.npy` both exist, but with nothing to place one in the other's frame the
  fan-out is not even measurable, and there is too little IF nucleus detection to fit one. 112 is also absent from `xenium_paths.csv` for a same-slide serial;
  only 6533's `0080870` is registered there.

Both are named rather than dropped. A donor whose registration failed supports **nothing** at
cell level, and this is the sentence that says so.

### 6539 is one stamp away

6539 has the best registration in the cohort (held-out p50 **0.524 µm**, the tightest of the
ten gated donors), a recorded channel-to-antibody mapping, and its IF values sitting on disk at
`data/development/restain/protein_rna/6539_0059865_*`. They were never written onto its Xenium
cell table, so there is presently nothing to join per cell. This is a gap in plumbing, not in
measurement, and it is the cheapest available expansion of the per-cell lane from two donors to
three.

Note also that **only 5 sections have a recorded channel-to-antibody mapping** at all
(`restain_sections.csv`: 6539, 6436, 6442, 115, 6591). That mapping exists in no file the
instrument wrote — the qptiff carries channel names but no biomarker tag. For the other 13
registered donors the IF channels cannot be *named*, so even a perfect registration would
produce values nobody can say the meaning of.

---

## 3. The many-to-one fan-out, quantified per donor

The IF image and Xenium segment independently. Measured nucleus-to-nucleus, through each
donor's **refined** transform, capped at 3.5 µm (half a nucleus diameter):

| donor | IF nuclei | Xenium nuclei | ratio | NN p50 (µm) | Xen nuc w/ IF partner | IF nuc w/ Xen partner | 1:1 share | fan-in max |
|---|---|---|---|---|---|---|---|---|
| 6436 | 395,048 | 696,464 | 0.567 | 1.534 | 59.1% | **98.3%** | 88.8% | 4 |
| 6442 | 584,723 | 958,918 | 0.610 | 1.278 | 63.5% | **99.0%** | 90.2% | 4 |
| 6539 | 410,469 | 291,026 | 1.410 | 0.675 | 74.0% | 51.1% | 94.7% | 3 |
| 6482 | 606,574 | 834,153 | 0.727 | 0.980 | 63.4% | 83.1% | 90.8% | 4 |
| 115 | 747,501 | 1,002,911 | 0.745 | 4.354 | 39.2% | 49.8% | 89.6% | 4 |
| 6493 | 326,092 | 658,912 | 0.495 | 4.474 | 46.6% | 89.8% | 90.6% | 3 |

(full table in `fanout.csv`; across the 16 measurable donors the 1:1 share spans
**81.4%–94.7%**, NN p50 **0.675–4.474 µm**, and fan-in never exceeds **4**.)

Three things follow, and the asymmetry is the point.

1. **The fan-out is modest and bounded.** Where a Xenium nucleus has an IF partner, that
   partnership is one-to-one for 81–95% of them and never more than four Xenium nuclei share one
   IF nucleus. This is not a segmentation catastrophe.
2. **The two segmentations disagree mostly by *count*, not by assignment.** The IF image finds
   0.50–1.41 nuclei per Xenium nucleus depending on donor. On 6436 and 6442 nearly every IF
   nucleus finds a Xenium partner (98.3%, 99.0%) while only 59–63% of Xenium nuclei find an IF
   one — the IF segmentation is the *sparser* of the two, and the missing pairs are Xenium
   nuclei the IF image did not resolve, not mismatches.
3. **None of this touches the IF values themselves.** They were measured in the Xenium cell's
   own polygon. The fan-out above is the ambiguity a *lifted call* would have to pass through,
   which is why section 4 reports the agreement both with and without it.

---

## 4. The direct per-cell protein-vs-RNA agreement

For **6436 and 6442 only**. Cells with no nucleus are excluded before anything is computed
(6436: 26,829 of 695,113; 6442: 51,392 of 958,744). Of the remainder, every cell fell inside the
imaged area — `n_outside_imaged_area` is 0 for both — and cells the RNA lane left `Unavailable`
are dropped, leaving **616,017** and **856,419** cells with both a protein value and an RNA call.

Agreement is reported as **rank AUROC**, threshold-free on purpose: naming a cut-off on
`feat__if_INS` would be asserting a positivity nobody has defined.

### Overall, and with the fan-out carried

| donor | stratum | n | INS AUROC (Endocrine) | GCG AUROC (Endocrine) |
|---|---|---|---|---|
| 6436 | all | 616,017 | **0.401** | **0.833** |
| 6436 | nucleus 1:1 | 347,368 | 0.419 | 0.839 |
| 6436 | shared IF nucleus | 33,508 | 0.429 | 0.777 |
| 6436 | unmatched | 235,141 | 0.383 | 0.837 |
| 6442 | all | 856,419 | **0.736** | **0.738** |
| 6442 | nucleus 1:1 | 524,628 | **0.778** | **0.773** |
| 6442 | shared IF nucleus | 45,043 | 0.620 | 0.624 |
| 6442 | unmatched | 286,748 | 0.687 | 0.703 |

**The fan-out costs about 0.15 AUROC, and it is visible.** On 6442 both hormones separate the
RNA-called Endocrine cells at 0.778 / 0.773 where the nucleus is matched one-to-one, and at
0.620 / 0.624 where it is shared. That gap is the price of segmentation ambiguity, measured on
the same cells rather than assumed. It is the number to carry onto anything derived from a
lifted call, and it says a many-to-one transfer is not free.

**On 6436 the gap is absent, because there is nothing there to degrade.** `feat__if_INS` sits at
or below chance in *every* stratum (0.383–0.429). A segmentation artefact would move with the
strata; this does not. The INS channel on section 0059907 carries no endocrine signal, while its
GCG channel is the strongest in the pair (0.833 overall; median 22.56 inside RNA-called Endocrine
against 1.41 outside, a 16-fold separation). `restain_sections.csv` already records Opal 570 as
FAILED on this section; this is an independent, per-cell measurement of a second channel
problem, and it is a statement about **this section's stain**, not about INS and not about the
method. You cannot indict a marker, only a method — and here not even that: one section.

### Within the endocrine call, on the types the RNA panel had to infer

INS, GCG and SST are **not in the 5K panel**. The RNA lane calls Beta/Alpha/Delta from surrogate
cores. The re-stain is the only place those inferred calls meet the hormone directly:

| donor | RNA fine type | n | INS AUROC vs other endocrine | GCG AUROC vs other endocrine |
|---|---|---|---|---|
| 6436 | Alpha | 14,346 | 0.330 | **0.805** |
| 6436 | Beta | 2,210 | 0.518 | 0.586 |
| 6436 | Delta | 982 | 0.506 | 0.403 |
| 6442 | Alpha | 4,408 | 0.580 | **0.816** |
| 6442 | Beta | 5,723 | **0.853** | 0.620 |
| 6442 | Delta | 822 | 0.609 | 0.509 |

GCG marks the RNA-called Alpha cells on both sections (0.805, 0.816) — two independent sections
agreeing, through different evidence, on the same claim. On 6442, where the INS channel works,
INS marks the RNA-called Beta cells at 0.853. On 6436, where it does not, the same comparison
returns 0.518, which is chance.

**Both directions are informative and neither side was tuned to improve the other.** No RNA call
was adjusted to raise an AUROC, and no IF value was rescaled. Delta separates on neither hormone
on either section, which is what the panel predicts — SST is absent from both the 5K panel and
the re-stain, so nothing here can speak to Delta and this table does not pretend to.

### What this calibrates, and what it does not

This is the project's **only** direct protein-vs-RNA agreement measurement, and it is the
reference the donor-level compositional comparisons should be read against. But its reach is
narrow and should be stated as such:

- **Two donors, two hormones, one lineage.** It calibrates the Endocrine axis. It says nothing
  about T_NK, B_Plasma, Myeloid or Granulocyte concordance, because no immune protein was
  re-stained. Post-Xenium tissue has been through protease digestion; INS and GCG are abundant
  secretory proteins and survive it, and lower-abundance surface markers such as CD3e reportedly
  do not.
- **No cross-donor quantity is computed here, so no eta-squared is reported.** With n=2 donors an
  eta-squared on `stain_group` or `acq_run` against a permutation null would be arithmetic, not a
  measurement, and reporting one would dress up a two-point contrast as a batch analysis. When
  6539 is stamped this becomes n=3, which is still not enough. The constraint is the design, not
  the code.
- **The agreement is between an IF value and an RNA-derived call, not between two calls.** The
  protein side of this comparison is a feature. It is not a PhenoCycler lineage call and it must
  not be quoted as one.

---

## 5. Constraints on anything built from this

**Cell identifiers do not cross modalities.** PhenoCycler `object_id` is a QuPath UUID; Xenium
`cell_id` is 10x segmentation, unique only within a bundle. There is no shared namespace. The
measured collision — 6442's 958,744 pancreas calls hitting 44 of its 177,799 pLN cells — is 44
id collisions, not a mapping, and joining on identifier would silently produce exactly that.
Everything above joins on **geometry through a stated transform**, never on an identifier.

**Do not co-plot frames.** PhenoCycler-to-Xenium rotations across the cohort span −168.9° to
+88.6°, with 12 of 20 donors above 20°. Four pixel scales are in play: Xenium 0.2125, IF 0.2505
(maintainer ruling), PhenoCycler 0.5077663810243286, DOMAIN 0.5078. No transform is applied in
the app; reprojecting 6442's `x_um` through the stored affine and differencing against
`x_um_reg` gives max 207 µm, mean 56.8 µm. Two layers drawn in one axis are two different
coordinate systems wearing the same units.

There is a related trap inside this lane. `warped_if_um.npy` is the **valis-stage** warp of
Xenium **cell centroids**, and it is wrong twice over for a nucleus comparison: a cell centroid
is displaced from its own nucleus, and the valis warp is superseded by the refined one. Measured
on 6442, the two warps differ by a mean of 1.37 µm and a maximum of **119 µm**. The fan-out in
section 3 warps Xenium *nuclei* with the *refined* transform for exactly this reason. An earlier
pass of this analysis used `warped_if_um.npy` and reported a median nearest-neighbour distance
of 3.56 µm for 6442 where the correct comparison gives 1.28 µm.

**Across tissues remains compositional only.** Nothing here reaches pancreas-to-pLN. A lymph
node has no islet, so a signed distance to one is undefined there — not zero. The maintainer's
forthcoming pLN-portions-of-pancreas-slides dataset may make an islet-to-lymphoid distance
meaningful for anatomically attached structures; that classification arrives as data and is not
inferred here.

---

## 6. Acceptance, restated as claims

1. **Per donor, what the linkage supports:** 2 per-cell (6436, 6442), 14 per-region, 2 nothing
   (6533 — registration failure at 311.5 µm median nucleus separation; 112 — no transform). Four
   per-region donors (6380, 6461, 6476, 6493) have a transform whose cell-level accuracy was
   never gated, and are per-region for that reason rather than on a measured failure.
2. **Fan-out quantified and carried:** 1:1 share 81.4–94.7%, fan-in never above 4, and the cost
   measured directly as **~0.15 AUROC** between one-to-one and shared-nucleus cells on 6442.
3. **Direct per-cell agreement where defensible:** 6442 INS 0.778 / GCG 0.773 on one-to-one
   nuclei; 6436 GCG 0.833 with INS at chance on a section whose stain is independently recorded
   as partly failed. Alpha is marked by GCG on both sections (0.805, 0.816); Beta by INS on 6442
   (0.853) and not on 6436 (0.518); Delta by neither, which the panel predicts.
4. **Where it is not defensible, named:** 6533 and 112 at cell level; all 14 per-region donors
   for any per-cell claim; every donor for any immune-lineage protein-RNA agreement, because no
   immune protein was re-stained.

Nothing above is a release. The typing latch has not moved.
