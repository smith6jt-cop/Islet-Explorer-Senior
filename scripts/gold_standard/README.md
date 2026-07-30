# Blinded broad-type audit

This directory provides a detached human audit of one exact, completed
Phenocycler typing run. It does not train or tune calibration, marker
thresholds, typing rules, donor exclusions, or release decisions.

The audit answers two separate questions:

1. How consistently can two blinded reviewers apply the broad-type
   vocabulary from the images?
2. On the cells where those reviewers agree, how do the frozen pipeline calls
   compare, including explicit `Ambiguous` and `Unavailable` abstentions?

Low reviewer agreement indicates that the annotation protocol is not a stable
reference. It is not permission to tune the pipeline against either reviewer.

## Vocabulary

Review each sampled cell as one of:

- `Immune`
- `Epithelial` — includes endocrine cells in this broad scheme
- `Endothelial`
- `Neural`
- `Mesenchymal`
- `Other` — a real cell that is none of the five broad types
- `Unsure` — the image does not support a confident human call

`Unsure` is preferable to a forced guess. The scorer reports it separately and
excludes it from pipeline accuracy, along with reviewer disagreements.

## 1. Build an immutable blinded bundle

First require a complete, current source run:

```bash
python -m phenocycler.pipeline status --config config.ini
```

Then snapshot a selected-donor audit bundle:

```bash
python scripts/gold_standard/build_bundle.py --config config.ini
```

The default output is:

```text
data/gold_standard_v2/<run_id>/
```

The builder validates `run.json` and the type-stage manifest, joins cells and
assignments by `(donor_id, object_id)`, and refuses to overwrite an existing
bundle. Predictions and their provenance are frozen in the private
`sample_key.csv`; scoring never reruns typing.

Sampling is stratified for audit efficiency:

| stratum | default share | definition |
|---|---:|---|
| `immune_called` | 40% | accepted frozen broad call is `Immune` |
| `islet_proximal` | 25% | remaining cell has QuPath-derived region `core` or `peri` |
| `uniform` | 35% | all remaining canonical cells |

The strata are disjoint. Each row has inverse sampling weight
`stratum_N / stratum_n`. Design-weighted results target the selected donors'
canonical-cell population only, not all donors in the cohort. Reviewer
disagreement and `Unsure` exclusions can limit that population interpretation
and are therefore reported alongside every score.

## 2. Annotate while blinded

Two reviewers must use separate QuPath projects, or separate copies of the
project data, and must not see:

- `sample_key.csv`;
- pipeline assignment overlays for sampled cells;
- the other reviewer's calls.

Run the importer on each donor image, passing the bundle's `qupath` directory:

```bash
QuPath script --project project.qpproj --image "exact image name" \
  scripts/gold_standard/import_gold_standard.groovy \
  --args "/absolute/path/to/data/gold_standard_v2/<run_id>/qupath"
```

In QuPath, sort detections by `gs_review_id`, inspect any channels needed, and
set each sampled cell to one vocabulary class. Export each reviewer's results,
passing the bundle directory and annotator ID:

```bash
QuPath script --project reviewer-A.qpproj --image "exact image name" \
  scripts/gold_standard/export_gold_standard.groovy \
  --args "/absolute/path/to/data/gold_standard_v2/<run_id>,A"
```

Use `B` for the second reviewer. The exporter writes separate files under
`completed/` and refuses to overwrite an existing image/annotator file.
Custom annotator IDs may contain only letters, numbers, `.`, `_`, and `-`.

## 3. Score the frozen calls

After both reviewers have exported every sampled image:

```bash
python scripts/gold_standard/score.py \
  --bundle data/gold_standard_v2/<run_id>
```

The scorer validates the bundle fingerprints and exact
review-ID-to-UUID mappings. It requires exactly one row per sampled cell from
each of the two declared reviewers, and refuses to score while either reviewer
has blank cells. Missing frozen predictions are fatal.

It reports:

- completion, `Unsure`, disagreement, raw agreement, and Cohen's kappa;
- weighted prediction coverage over the full frozen sample;
- consensus-conditional weighted accuracy with abstentions counted as uncalled;
- consensus-conditional weighted accuracy when the pipeline makes a call;
- weighted broad-type precision, recall, and F1;
- per-donor coverage and accuracy;
- raw and design-weighted confusion matrices retaining `Ambiguous` and
  `Unavailable`.

The first successful score creates `scores/` and will not overwrite it.
`scores/score_manifest.json` fingerprints every completed annotation input and
score output and records the bundle and source-run identities.

## Bundle contents

| file | purpose | reviewer access |
|---|---|---|
| `manifest.json` | immutable run, typing, sampling, and file provenance | yes |
| `sample_key.csv` | UUID, stratum, weight, and frozen pipeline prediction | **no** |
| `blinded_review_index.csv` | read-only review ID and location index | yes |
| `qupath/gold_standard_<donor>.csv` | UUID-to-review-ID importer input | yes |
| `completed/annotations_<reviewer>_<image>.csv` | exported human calls | respective reviewer |
| `scores/` | immutable, fingerprinted audit summaries | only after both reviews finish |

Existing `data/gold_standard_v1` material is historical and must not be
overwritten. If external v1 annotation has begun, preserve it and create a v2
bundle for the current canonical cell universe.
