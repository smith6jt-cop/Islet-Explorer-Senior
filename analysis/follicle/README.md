# Spleen follicle analysis

This is the consolidated SpleenFollicleCounterQP lane. It contains the tracked
QuPath workflow, region-layer analysis, 20 hypothesis notebooks, and paired
H&E/DAPI follicle U-Nets. Large images, QuPath objects, measurement exports,
and generated results remain outside Git.

## Layout

```text
analysis/follicle/
├── analysis/                 notebooks and shared data_utils.py
├── classifiers/              reviewed QuPath classifiers
├── config/sample_manifest.csv
├── resources/                QuPath displays
├── scripts/                  producers, validation, and U-Net workflows
├── tests/                    synthetic geometry and U-Net tests
├── identity.py               exact alias-to-sample resolution
├── paths.py                  code/data/output root contract
├── project.qpproj            tracked QuPath descriptor
└── IMPORT_PROVENANCE.md      source-repository evidence
```

Runtime roots are explicit:

- Code: this directory.
- Data: `FOLLICLE_DATA_ROOT`, default `/data/SpleenFollicleCounterQP`.
- Output: `FOLLICLE_OUTPUT_ROOT`, default the data root.

The manifest is authoritative for real sample identity. Unknown aliases abort;
the workflow never derives donor identity from a filename. Generated U-Net
metadata carries an explicit `sample_id` for the same reason.

## Setup and checks

From the monorepo root:

```bash
python -m venv .venv
.venv/bin/pip install -r analysis/follicle/requirements.txt

.venv/bin/python analysis/follicle/scripts/preflight.py --profile all
MPLCONFIGDIR=/tmp/islet-mpl .venv/bin/python -m pytest -q analysis/follicle/tests
MPLCONFIGDIR=/tmp/islet-mpl .venv/bin/python \
  analysis/follicle/scripts/verify_analysis_state.py --static-only
```

The preflight reads no large files. It lists the complete missing set and exits
2 until the external contract is restored deliberately. At the 2026-09-03
import checkpoint, `Images/` and QuPath `data/` pass, while the four required
files under `Measurements/` are absent.

## Main workflow

Once preflight passes, run producers in dependency order:

```bash
.venv/bin/python analysis/follicle/scripts/fit_area_floors.py
.venv/bin/python analysis/follicle/scripts/radial_follicle_profiles.py
.venv/bin/python analysis/follicle/scripts/build_region_layer_qupath.py
.venv/bin/python analysis/follicle/scripts/validate_region_layer.py
.venv/bin/python analysis/follicle/scripts/compute_primary_endpoint.py
.venv/bin/python analysis/follicle/scripts/plot_region_layer_results.py
analysis/follicle/scripts/run_notebooks.sh
MPLCONFIGDIR=/tmp/islet-mpl .venv/bin/python \
  analysis/follicle/scripts/verify_analysis_state.py
```

Outputs land below `$FOLLICLE_OUTPUT_ROOT/analysis/`. Notebook results treat
donor as the independent unit; cells, follicles, images, and sampled windows
remain nested observations. See [docs/analyses.md](docs/analyses.md) for the
hypothesis inventory and [docs/data.md](docs/data.md) for the legacy data
details. Where those historical documents mention repository-relative data,
the explicit root contract above takes precedence.

## U-Net workflow

`scripts/export_HE_training_pairs.groovy` and
`scripts/export_DAPI_training_pairs.groovy` export aligned images, masks, and
metadata. Python trains, predicts, and reports:

```bash
.venv/bin/python analysis/follicle/scripts/train_HE_unet.py --folds 4
.venv/bin/python analysis/follicle/scripts/train_DAPI_unet.py --folds 4
.venv/bin/python analysis/follicle/scripts/check_dapi_parity.py --strict
```

Root workflow `.github/workflows/follicle-ci.yml` exercises both modalities on
synthetic donors, including pooled training and cross-modality warm start. The
smoke scores are crash checks only; they are not biological performance claims.

H16/H20 outputs are retained as historical evidence. They do not declare the
Phenocycler identity producer for the unified islet workflow.
