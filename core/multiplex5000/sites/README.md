# Site files

A site file says where a machine keeps its **data**, so that which checkout the code is installed
from never decides which data it reads. `multiplex5000/roots.py` reads it; the
`phenocycler` import remains a compatibility alias.

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

The [HiPerGator Xenium source index](hipergator/README.md) records exact paths
for 65 Xenium region outputs. A configured HiPerGator site names it as
`sources.xenium_source_index`; Stage 04 binds it and each selected raw matrix
and cell table in its receipt. A local pilot may instead use the accepted
Stage 00 paths and `sources.xenium_records/MANIFEST.tsv`.
`settings.rna_acquisitions` is a legacy configuration key selecting explicit
inventory IDs; those rows are region outputs, not instrument-run counts. Keep
run, slide, section and region identifiers distinct. A single run or slide can
contain several regions.

The dated full-cohort HiPerGator release is under
`/blue/maigan/smith6jt/mplx5k-io60-rna-20261027/releases/io60-candidate-1`;
`site-current.json` and `venv-current/bin/python` under that project are recorded
runtime locations. See the [current evaluation guide](../../../docs/PORTABILITY_REMEDIATION_USER_GUIDE.md)
for the evidence and scope. This corrects the older statement that HiPerGator
had no Stage 00 release: do not repeat its setup or RNA execution on that basis.
Before a new Stage 04 execution, its own selected release still needs the
required accepted Stage 00 output. The two-region workstation pilot and
65-region HiPerGator execution remain separate processing records. Recorded
remote paths do not establish deployment of today's workstation source; verify
that site's actual producer/runtime before using it.

The [dedicated pLN v3 export index](pln_qupath_exports_v3.csv) binds each of the
19 recorded images in `pln_node_annotations.csv` to one QuPath cell-export filename
and segmentation ID. Its filenames are relative to `sources.pln_v3_export_cells`,
which each site sets to the folder holding those Parquet files. The source loader
checks the roster identity, QuPath schema, first `Image` value, and export row
count. The roster's older counts are comparison evidence from a prior
segmentation; a mismatch does not reject a later export. This binding is a
source check, not compensated protein measurement or typed cells.
Run `mplx5k sources pln` on a configured site, or
`mplx5k sources pln --root /absolute/path/to/export/cells` to inspect an
explicit folder. `--image` limits it to recorded images for a bounded pilot.
This command reads Parquet metadata and the first `Image` value; it does not
claim every cell in a large file was checked or accept a processing gate.

## The MultiPlex5000 environment

**Current workstation use (November 6, 2026):**
`/data/envs/multiplex5000/bin/python` already exists and imports live Senior
`multiplex5000` plus its `phenocycler` compatibility alias from a neutral working
directory. Its `mplx5k` script may not be on the shell's PATH; selecting the
interpreter explicitly avoids that dependency:

```bash
export MPLX5K_PYTHON=/data/envs/multiplex5000/bin/python
export MPLX5K_SITE_FILE=/home/smith6jt/.config/multiplex5000/site.json
env -u PYTHONPATH "$MPLX5K_PYTHON" -m multiplex5000.cli context \
  --project io60_xenium_pilot --releases-root /data/io60-current/releases
```

`context` reports the resolved interpreter, package, producer/site/project and
output bindings without executing a stage. It does not rehash all scientific
inputs or establish biological acceptance. Preserve the reported context when
recording a new analysis. Explicit project and release-root arguments belong
on subsequent run, status and gate commands as well.

The repository `.venv` remains separate: the same neutral-directory check found
no `multiplex5000` there and `phenocycler` from Senior-v3. A `PYTHONPATH` override
can change which code that interpreter imports, but cannot make it the same
runtime as `/data/envs/multiplex5000`. Keep existing historical environments and
their consumers intact. Pinned CODEX replay and HiPerGator/Arrow runtimes remain
separate contracts; this workstation choice does not replace them.

### Historical environment creation recipe

The following D1 recipe explains how the workstation environment was created.
It is not a command to rebuild the existing environment. For a new installation,
choose a new path, review the dependency requirements for the selected workflow
and verify the resulting runtime before processing data.

```bash
# Historical creation at a then-new path; do not rerun over an active environment.
/usr/bin/python3.12 -m venv /data/envs/multiplex5000
/data/envs/multiplex5000/bin/python -m pip install "pip==26.2.1"
/data/envs/multiplex5000/bin/python -m pip install -r requirements-multiplex5000.lock
/data/envs/multiplex5000/bin/python -m pip install --no-deps -e core/multiplex5000
make env-check PY=/data/envs/multiplex5000/bin/python
```

`requirements-multiplex5000.lock` is `requirements.lock` without the three pins
its header names. The retired dct lane and the viewer's historical runtime are
not a reason to alter the checked module environment. For a compatible current
Make target, pass the interpreter explicitly, for example
`make test-fast PY=/data/envs/multiplex5000/bin/python` from the repository root.
