# DeepCellTypes pilot — per-image review protocol

```yaml
status: SUPERSEDED
superseded_by: DEEPCELL_TYPES_REVIEW_2026-08-24.md
superseded_on: 2026-08-24
superseded_because: >
  Its two load-bearing assumptions were measured and found false -- max-softmax does not
  predict agreement, and the reviewed object (the argmax call) discards what the probability
  vector carries. No review field was ever opened under it and no reviewer was ever named,
  so nothing it governed was spent. The body below is UNEDITED on purpose: this document
  said a change of mind must be a new dated protocol, and its integrity is the point.
written: 2026-08-21
creates_cell_calls: true      # the first data/development artifact in this repo that does
owner: smith6jt
reviewer: TBD — record the name here before the first field is opened
lane: mutable, data/development/. Writes no StageManifest and does not mutate data/runs/.
```

**Written before any inference has run.** Nothing below may be edited after the first review
field is opened; a change of mind is a new dated protocol that supersedes this one and says so.

Every other review artifact in this repo carries `creates_cell_calls: false` — the fundamental-pair
lane is explicit that "it does not call markers or create cell types." This one does. That is the
whole reason it needs a bar fixed in advance.

## What is being decided

Whether DeepCellTypes' per-cell calls on this cohort are right enough to build on. **This
per-image verdict is the release gate.** Nothing downstream of it happens in this piece of work:
no cohort run, no `apps/islet` change, no promotion to a `vision_types` stage.

A negative verdict is a real, expected outcome, and the "How to discard this" section below
exists so that acting on it is a single command rather than an archaeology project.

## Field selection — fixed rule, not scrolling

Per donor: **12 fields of 500 × 500 px** at level 0, selected by a rule that cannot see the
calls' plausibility.

1. Enumerate every 500-px window containing **≥ 40 cells** with a call.
2. Compute per-donor max-softmax **terciles** over all called cells. Record the two boundary
   values in this file's results appendix *before* reviewing. These are the three **confidence
   strata**: low / middle / high.
3. Assign each window the stratum of its own median max-softmax.
4. Within each stratum, order candidate windows by `sha256(f"{donor}:{y0}:{x0}")` ascending and
   take the first **3**. That is deterministic, reproducible from the donor ID alone, and blind
   to content.
5. Add **3 islet-anchored fields**: the same hash ordering, restricted to windows containing at
   least one cell whose QuPath parent is `Islet_N`. The study question is the peri-islet immune
   environment; without this, a hash draw over whole-pancreas area can miss islets entirely.

Record the 12 window origins per donor. If a window cannot be rendered, take the next one in the
same hash order and record the substitution and why.

## What the reviewer records

For every cell they adjudicate, one of:

| Verdict | Meaning |
|---|---|
| `Correct` | The specific class is right |
| `Broad-only` | Wrong specific class, right broad compartment (e.g. CD4T called CD8T) |
| `Wrong` | Wrong broad compartment |
| `Cannot tell` | The image does not carry enough evidence either way |

`Cannot tell` is a first-class answer, not a failure to decide. It is excluded from the
denominator and reported separately — a stratum where most cells are `Cannot tell` has not been
reviewed, whatever its pass rate says.

### Per-class tally sheet

One row per specific class actually called in the reviewed fields, per donor, per stratum:
class, cells reviewed, Correct, Broad-only, Wrong, Cannot tell. Blank rows for classes that never
appeared are kept — a class that was never called is a result.

## The bar, stipulated in advance

Fixed here on 2026-08-21, before any output exists. Its one defensible property is that it was
chosen before anyone had seen a result.

**PASS for a donor requires, in the high-confidence stratum:**

- `Wrong` ≤ **10 %** of adjudicable cells (Correct + Broad-only + Wrong); and
- no single specific class accounts for > **25 %** of all `Wrong` calls — a concentrated error is
  a systematic one and is not averaged away; and
- ≥ **50 %** of reviewed cells adjudicable (not `Cannot tell`), or the stratum is reported as
  not reviewed.

**The verdict is recorded separately for each of the three strata.** A pass in the high stratum
does not license the middle or low ones; it bounds where the calls may be used to exactly the
stratum that earned it. There is no aggregate pass across strata and none may be computed later.

**Lane PASS requires all three pilot donors to pass.** Two of three is a fail with a finding, not
a partial pass.

These are acceptance bars on human adjudication. They are not marker thresholds, and nothing in
this protocol may be tuned toward an expected cell-type frequency, cross-donor uniformity, or
disease ordering.

## Blinding — partial at best, stated plainly

The maintainer stained and scanned these slides and knows the donor IDs, the disease labels and
the batch structure. Full blinding is not available and pretending otherwise would be worse than
naming the limit. What **is** achievable and is required:

- The reviewer sees the specific call and the confidence stratum, **not** the model's
  probability vector and **not** which stratum is "supposed" to be good.
- Field order is shuffled across strata and donors so the strata are not reviewed in blocks.
- Disease status is not displayed in the review bundle.
- The three confidence-stratum boundary values are recorded before review and not revisited.

Record any further un-blinding that occurs, rather than discarding the field.

## Opening state of the bundle

QuPath's `view_phenotypes.groovy` starts with **every** class visible — with 51 specific classes
plus `Ambiguous` and `Unavailable` that is a ~53-item checkbox list over a solid carpet of dots,
and any verdict formed from it would be a verdict about the densest classes, not the immune ones
the study asks about. So the bundle ships with a defined opening state:

- **broad level first**, immune classes on, everything else hidden;
- a second CSV restricted to the confidence-stratified, spatially-sampled subset above, so the
  reviewer's fields come from the rule and not from scrolling.

## How to discard this cleanly

If the verdict is negative, everything this lane produced is removed by deleting these paths.
Each is listed so a negative verdict is executable rather than aspirational. **Nothing outside
this list is written**, and in particular nothing under `core/phenocycler/data/runs/`.

```text
data/development/deepcell_types_review/          # per-donor outputs, figures, review bundle
core/phenocycler/scripts/dct/                    # pilot code (git revert, do not rm -rf)
core/phenocycler/docs/DEEPCELL_TYPES_REVIEW_*.md # this protocol and its results appendix
core/phenocycler/pyproject.toml                  # remove the [dct] optional-dependency group
~/.deepcell/models/                              # cached checkpoints, if no longer wanted
```

`git reset --hard` and `git clean -xfd` are **unsafe in this repository** — 157.9 GB of untracked,
ignored data lives under `core/phenocycler/data/`. Roll code back with `git revert`.

## Results appendix

Left empty on purpose. Filled in only after review, in this file, below this line.

---
