# Site files

A site file says where a machine keeps its **data**, so that which checkout the code is installed
from never decides which data it reads. `phenocycler/roots.py` reads it (MultiPlex5000 plan, phase
D1).

- **Where it is looked for:** the file named by `$MPLX5K_SITE_FILE`, else
  `~/.config/multiplex5000/site.json`. With neither, the data roots fall back to the checkout's own
  `data/` and a `SiteFallbackWarning` says so — the right behaviour for CI and a fresh clone.
- **What is refused:** a named file that is missing, unreadable JSON, a different `schema`, a
  relative `data_root`, or a `data_root` that is not a directory. An explicit setting that cannot
  be honoured is an error, never a fallback.
- **Format** (`"schema": "multiplex5000-site/1"`):

| key | required | meaning |
|---|---|---|
| `site` | yes | the machine's name, reported by `make env-check` |
| `data_root` | yes | absolute path of the repository-level `data/` tree (`data/development` lives under it) |
| `layers_dir` | no | absolute path of the app layer artifacts; `$ISLET_LAYERS_DIR` still wins when set |
| `sources` | no | recorded sources that live outside the repository, by name: `donor_notes` (the maintainer's donor notes; only their donor and post-Xenium IF sections are read), `donor_notes_copy` (the repo-root copy, compared with it), `pln_blocks` (`results/v3/pln_blocks.csv`), `qupath_manifest` (the generated, untracked QuPath manifest in the main checkout), `pancreas_qptiffs` and `pln_qptiffs` (the folders of both PhenoCycler projects' qptiffs; stage 00 reads each one's first page for its panel, experiment and scan time) and `xenium_records` (copies of the `experiment.xenium`, `gene_panel.json` and `metrics_summary.csv` of Xenium outputs this machine holds only in part, taken from the complete outputs on HiPerGator, with `MANIFEST.tsv` giving each file's source and sha256), `donor_metadata` (the project's `donor_metadata_panc.xlsx`; stage 00 reads each donor's `disease.status`) and `npod_cases` (nPOD's case export; its `donor_type` gives the status of a donor the metadata lacks, and is checked against it where both have the donor). Read by `multiplex5000.relationships` and stage 00; a source the site does not name is unknown, not absent |
| `releases_root` | no | where releases are built, one folder per release (D3); default `/data/io60-current/releases`; `$MPLX5K_RELEASES_ROOT` still wins when set |
| `snakemake` | no | the `snakemake` of the graph's own environment (`requirements-graph.txt`), for `mplx5k graph`; else the one on PATH |
| `tmp_dir` | no | where every stage run puts its temporary files (`TMPDIR`, which pool workers inherit) and DuckDB's overflow, in a folder of the run's own that is removed when the run ends; the workstation's is on `/data`, because its home disk is small. Else Python's own temporary folder |
| `note` | no | free text |

Files here are templates, one per machine; copy the right one to `~/.config/multiplex5000/site.json`.
They hold paths only — no accounts, credentials or scientific settings. Tracked inventories
(`core/phenocycler/data/integration`) are not site data: they are versioned with the code and follow
the checkout.

| file | machine |
|---|---|
| `sn4622120742.json` | the workstation |

A HiPerGator file is written when the HPG lane moves onto MultiPlex5000 (phase F).
The separate [HiPerGator Xenium source index](hipergator/README.md) records
exact paths for the 65 accepted Xenium region outputs. A future HPG site file
must name it as `sources.xenium_source_index`; Stage 04 binds it and each raw
matrix and cell table in its receipt. The accepted Stage 00 directory must
also be present under that project's release root before Stage 04 runs.
It is not currently on HiPerGator, so do not treat a local Stage 04 test as
a full-cohort HPG ingest.

## The MultiPlex5000 environment

The repo `.venv` is shared by several checkouts and, since 2026-10-03, imports `phenocycler` from
the v3 worktree (`make env-check` reports WRONG CHECKOUT). Rather than reinstall it under the
sessions still using it, phase D1 builds a second environment at a new path, with this checkout
installed editable:

```bash
/usr/bin/python3.12 -m venv /data/envs/multiplex5000
/data/envs/multiplex5000/bin/python -m pip install "pip==26.2.1"
/data/envs/multiplex5000/bin/python -m pip install -r requirements-multiplex5000.lock
/data/envs/multiplex5000/bin/python -m pip install --no-deps -e core/multiplex5000
make env-check PY=/data/envs/multiplex5000/bin/python
```

`requirements-multiplex5000.lock` is `requirements.lock` without the three pins its header names
(the retired dct lane, which cannot be installed from it, and the TissUUmaps viewer fork, which
stays in the repo `.venv`). Run anything in this repository with this interpreter by passing
`PY=` to `make`, for example `make test-fast PY=/data/envs/multiplex5000/bin/python`.
