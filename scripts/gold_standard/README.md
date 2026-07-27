# Gold-standard annotation protocol (v1)

## Why this exists

There is no accuracy measurement anywhere in this project. Every QC in the pipeline is an
*internal-consistency* check — pair states, seed stability, divisor reproducibility, the compartment
accounting — and internal consistency is precisely what a systematic error survives. A ~55 pp
epithelial error passed 217 tests, 11 review bundles and three gates without one of them firing.

Nothing in the repo can currently answer **"what fraction of the cells called Immune are immune
cells?"**, and that is the question the study's conclusion rests on.

This bundle is the instrument for answering it. It measures nothing on its own.

## What you are being asked to do

Assign one of seven labels to each sampled cell, from the image alone:

| label | what it means here |
|---|---|
| `Immune` | any leukocyte — T, B, plasma, macrophage, monocyte, DC, granulocyte, mast, NK |
| `Epithelial` | acinar, ductal, **and endocrine** (islet cells are epithelial in this scheme) |
| `Endothelial` | vessel lining |
| `Neural` | nerve fibre, ganglion, Schwann |
| `Mesenchymal` | fibroblast, stellate, smooth muscle, pericyte, adipocyte |
| `Other` | a real cell that is genuinely none of the above |
| `Unsure` | you cannot tell from the image |

**`Unsure` is a legitimate answer.** A forced guess is worse than a recorded uncertainty: it
manufactures agreement where there is none and inflates the accuracy estimate. The scorer reports the
unsure rate separately and excludes those cells rather than adjudicating them.

Two rules that decide most edge cases:

- **Endocrine is Epithelial.** Islet cells go in `Epithelial`, not `Other`, even though they look
  nothing like acinar tissue. This is the pipeline's scheme and the comparison is meaningless if the
  two use different definitions.
- **When a cell could be two things, pick what it *is*, not what it is near.** A leukocyte sitting
  inside an islet is `Immune`. A capillary threading through an islet is `Endothelial`.

## How to do it

1. Run `import_gold_standard.groovy` on the donor's image (edit `folder` first). Sampled cells get
   PathClass `GS_review` and a `gs_review_id` measurement.
2. **Analyze → Show detection measurements**, sort by `gs_review_id`, and work down the list. Use
   whatever channels you need; there is no restriction on what you may look at.
3. Set each cell's class from the vocabulary above.
4. Run `export_gold_standard.groovy` (edit `annotator` to `A` or `B`) → writes to
   `data/gold_standard_v1/completed/`.
5. When both annotators are done: `python scripts/gold_standard/score.py --arm <arm> --design union`.

## Blinding — do not break these

- **The two annotators must not see each other's calls.** Work in separate QuPath projects or separate
  copies of the `.qpdata`. Inter-annotator κ is the number that bounds every accuracy claim, and it is
  worthless if the second annotator has seen the first one's answers.
- **Neither annotator may look at the pipeline's calls** for these cells before annotating. The
  annotation sheets and the QuPath import deliberately carry location only — no compartment, no
  stratum, no marker value. Do not run `import_broad_lineage.groovy` on the same project.
- Rows are shuffled, so donor and stratum cannot be read off the ordering.

## How the sample was drawn (and why it is not uniform)

Five donors spanning all three disease groups, all three acquisition-software versions, and both
batch-2 panel variants — so a per-compartment accuracy estimate is never confounded with the
ND ↔ Fusion-2.4.1 aliasing (Fisher p = 0.0147; see `data/donor_acquisition_batch.csv`).

500 cells per donor in three disjoint strata:

| stratum | share | why |
|---|---|---|
| `immune_called` | 40% | immune carries the study's conclusion and is a few percent of cells; a uniform sample would contain too few to estimate immune **precision**. Defined as CD3e/CD68/CD11b positive under *any* divisor arm — the widest net, so every candidate the pipeline could produce is eligible for review. |
| `islet_proximal` | 25% | `dist_islet ≤ 20 µm`. The islet boundary is where the study's signal lives and where the compartments are hardest to separate. |
| `uniform` | 35% | simple random over everything else. **This is what makes recall estimable** — without it the sample contains only cells the pipeline already called, and a population the pipeline misses entirely would score perfectly. |

Each cell carries `weight = stratum_N / stratum_n`, so the scorer reports both the raw sample rates and
population-weighted rates. **Quote the weighted ones**; the raw ones describe the sample, which was
deliberately built not to look like the cohort.

## Reading the result

The scorer prints three blocks in order of authority:

1. **Inter-annotator agreement.** If Cohen's κ < 0.6 the label definitions disagree, not just the
   calls — fix this protocol and re-annotate. Do not quote accuracy numbers computed under a low κ.
2. **Per-compartment precision / recall / F1** against the consensus (both annotators agree, neither
   said `Unsure`). Disagreements are excluded and counted, never adjudicated by majority — with two
   annotators that is just picking one — and never adjudicated against the pipeline, which would make
   the whole measurement circular.
3. **The confusion matrix**, which is where the errors actually go. An Epithelial cell called `Other`
   and an Epithelial cell called `Immune` are very different failures for this study.

## Files

| file | contents | give to annotators? |
|---|---|---|
| `sample_key.csv` | review_id, donor, object_id, coordinates, stratum, weight, pipeline flags | **no** |
| `annotation_sheet_A.csv` / `_B.csv` | review_id + location + blank label columns | yes (or use QuPath) |
| `qupath/gold_standard_<donor>.csv` | object_id → review_id for the importer | yes |
| `manifest.json` | seed, donors, strata, vocabulary, weighting, blinding rules | yes |
| `completed/annotations_<who>_<image>.csv` | written by the export groovy | — |
| `scores/` | written by `score.py` | — |
