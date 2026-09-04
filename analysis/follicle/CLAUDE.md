# Spleen follicle lane

This directory is the tracked SpleenFollicleCounterQP implementation inside
Islet Explorer Senior. It is an analysis lane, not a second repository or a
source of cross-modal cell-identity truth.

## Roots and entry points

- `CODE_ROOT` is this directory. It owns Python/Groovy source, notebooks,
  classifiers, the QuPath project descriptor, small reference tables, and
  `config/sample_manifest.csv`.
- `DATA_ROOT` is external and defaults to `/data/SpleenFollicleCounterQP`.
  Override it only with `FOLLICLE_DATA_ROOT`.
- `OUTPUT_ROOT` defaults to `DATA_ROOT`. Override it with
  `FOLLICLE_OUTPUT_ROOT`; generated regions, figures, tables, models, and
  predictions belong there, not in the tracked code tree.
- Run `python analysis/follicle/scripts/preflight.py` from the monorepo root
  before a real-data run. It reports every missing path together and exits 2
  on an incomplete contract.

As absorbed on 2026-09-03, the external `Images/` and QuPath `data/` trees are
present, but `Measurements/` is absent. Four required exports are currently
missing at their declared locations: `AnnotationsFinal.csv`, `ForSH2B3.csv`,
`Cells.csv`, and `Spleen_Upload2_2024Sept17_v2_distances.csv`. Do not restore,
move, or regenerate them implicitly.

## Identity and analysis invariants

1. `config/sample_manifest.csv` is the authority for every real image,
   measurement, region stem, and metadata alias. `identity.resolve_sample_id`
   performs an exact lookup and fails on unknown or duplicate aliases.
   Filename regexes, prefix extraction, and ALT-ID digit slicing are forbidden.
2. Synthetic U-Net datasets declare `sample_id` in each metadata JSON. Cross-
   validation groups on that field; multiple images from one donor must carry
   the same value.
3. Donor is the independent statistical unit. Image, follicle, cell, and
   window observations are nested within donor; larger images increase
   precision and are not extra donors.
4. Missing measurements are unavailable, never zero. `HDL018`, `HDL021`, and
   `HDL172` remain excluded unless a documented analysis explicitly includes
   the first two for its historical robustness check.
5. Region membership is exact point-in-polygon against the generated region
   layer. The QuPath `Parent` string is used only by the independent validation
   path. Pixel sizes are centralized in `analysis/data_utils.py`.
6. H16 and H20 are preserved historical evidence for the protein-producer
   comparison. Their RESTORE/CD20 labels are not a selected or released
   Phenocycler identity method and must not populate the unified cell contract.
7. Keep the H&E and DAPI U-Net scientific transforms paired with their Groovy
   exporters. Any transform edit requires parity tests. A smoke run proves
   execution, not model quality.

## Validation

```bash
MPLCONFIGDIR=/tmp/islet-mpl .venv/bin/python -m pytest -q analysis/follicle/tests
MPLCONFIGDIR=/tmp/islet-mpl .venv/bin/python \
  analysis/follicle/scripts/verify_analysis_state.py --static-only
.venv/bin/python analysis/follicle/scripts/preflight.py --profile all
```

The first two commands must pass without external data. The last command must
pass before notebooks or real-data producers run. Root CI additionally trains
both U-Nets on synthetic data, tests pooled and warm-start modes, verifies
byte-exact DAPI parity, predicts slides, and builds the concordance report.

## Import provenance

The exact source history and preservation branches are recorded in
`IMPORT_PROVENANCE.md`. The imported source tree was commit `1ed38d2`; its
upstream base was `6ae30ba`, and the pre-consolidation local-work archive was
`f6cbcd4`. Do not edit or delete the source repository as part of normal work
here.
