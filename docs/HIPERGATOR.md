# Running this on HiPerGator — setup and the basic analyses

A working outline for the 20-donor cohort: PhenoCycler protein and Xenium 5K RNA, pancreas
and pancreatic lymph node, on UF Research Computing.

Commands are the real ones (checked against `argparse`); paths are the real ones (checked
against `PipelineConfig`). Where a number is yours to supply rather than the repo's, it says so.

> **RESTORE is the production marker-calling method.** The content-addressed
> workflow retains its immutable `calibrate` stage key, but that stage now runs
> the validated RESTORE v10 pair fit and emits tri-state marker evidence. The
> intervening donor-calibration replacement was not an intended migration.

---

## The shape of it

```
OFF-CLUSTER                       ON HIPERGATOR  (one job, eight stages)
────────────                      ─────────────────────────────────────
QuPath segmentation               ingest      QuPath CSV  -> cells/
  ├─ CellMeasurements.csv ──┐     geometry    segmentation QC
  ├─ cells__<image>.geojson ├───▶ redsea      lateral spillover compensation
  └─ <image>.qptiff ────────┘     expression  measurable-marker selection
                                  controls    retained reference-control diagnostics
Xenium runs ──────────┐           calibrate   RESTORE pair fit       -> marker_evidence/
  └─ output bundles   │           type        ordered residual typing -> assignments/
                      │           states      process states         -> cell_states/
                      │                            │
                      └──────────▶ integration  ◀──┘   (second conda env, see Phase 5)
```

Everything upstream of `ingest` is declared **once** in a cohort manifest
(`data/qupath_manifest.json`). There is no per-stage script chain, no array range and no
hand-maintained donor count: the manifest is the donor set, and stage manifests carry
restart/current-state logic.

Two units matter and are easy to confuse:

- **Section** — one image, one GeoJSON, one coordinate frame. `6539` (pancreas) and `6539pln`
  (lymph node) are *two* sections. `cfg.discover_sections()` returns these.
- **Donor** — `6539`. `cfg.discover_donors()` parses section keys to donor ids, dedupes, and
  drops the cohort exclusions. This is what joins to metadata and to the Xenium manifest.

**Cohort exclusions are fail-closed** (`phenocycler/cohort.py`). Two donors are hard-excluded
and `ensure_eligible_donors` *raises* if either appears:

| donor | reason |
|---|---|
| `6457` | poor DAPI staining prevented accurate PhenoCycler channel registration |
| `6579` | pancreatitis outlier; excluded by maintainer decision |

That is not a filter you can pass around — it is enforced at discovery and again at stage
entry, so a run cannot quietly include them.

> ⚠️ **A known cohort constraint.** `artifacts._validate_cohort_images` requires **exactly one
> image per donor**. A donor with both a pancreas and a lymph-node section is therefore not
> yet representable in the manifest. Until `QuPathImageManifest` grows a section/ROI field,
> the pancreas and pLN halves must be built as **separate cohort manifests and separate
> runs**. Doing that field properly changes `content_id` and hence every `run_id`, so land it
> before a production run rather than after.

---

## Phase 0 — what happens off-cluster

Neither of these can be done from a SLURM job.

**0.1 QuPath, per image.** Segment, then export:
- the per-cell measurement table (`CellMeasurements.csv`);
- full-resolution boundaries → `cells__<image>.geojson`, via
  `scripts/groovy/export_cells_geojson.groovy`;
- optionally nucleus boundaries, if you want the nuclear/cytoplasmic geometry QC.

Edit the Groovy script's output folder first — it hardcodes `~/Phenocycler_Analysis/data/...`,
which is not where your data will live on HiPerGator.

You will need, per image, the **exact** QuPath `Image` column value and the pixel size. Both
go in the manifest verbatim; the manifest fingerprints every input file, so a renamed or
re-exported file is detected rather than silently used.

**0.2 Xenium bundles.** You need the output bundle per run (`cells.parquet`,
`cell_boundaries.parquet`, transcripts) — the polygons ship there and *not* in the zarr or a
slimmed h5ad. Note where any phenotyped h5ads from `Xenium_Analysis` live; Phase 5 needs that
`data/processed` root.

---

## Phase 1 — one-time setup

```bash
git clone --recurse-submodules https://github.com/smith6jt-cop/Phenocycler_Analysis.git
cd Phenocycler_Analysis
bash setup.sh                    # submodules + conda env + editable install + verification
conda activate phenocycler-analysis
```

`setup.sh` initialises submodules, creates/updates the env from `environment.yml`, runs
`pip install -e '.[test]'`, and then verifies that the packaged marker registry,
ordered residual rules, and optional probabilistic rules are installed and
that the registry/rule contract loads. If that verification fails, stop
— nothing downstream will work.

RESTORE no longer requires the historical external submodule or `spams`; the
validated v10 implementation is self-contained in `restore_validation.py` and
uses scikit-learn NMF. `external/XeniumPanelExplorer` remains the vendored
dependency used by the integration crosswalk.

**Put the data somewhere with room.** `config.ini` ships `data_dir = data`, i.e. inside the
repo, which on HiPerGator means your home quota:

```bash
export PHENOCYCLER_DATA_DIR=/blue/<group>/<user>/phenocycler/data
export PHENOCYCLER_QUPATH_MANIFEST=/blue/<group>/<user>/phenocycler/qupath_manifest.json
export PHENOCYCLER_RESTORE_COMPARTMENT_SAMPLE=/blue/<group>/<user>/phenocycler/qupath_compartment_full_20donor
export TMPDIR=/blue/<group>/<user>/tmp      # DuckDB spills here during ingest
```

Conda envs default to `~/.conda/envs`, also on the home quota; `conda config --add envs_dirs
/blue/<group>/<user>/envs` first if that is tight.

Budget: the QuPath CSV is large (~138 GB for the full cohort), the parquet output is tens of
GB, and REDSEA scratch (masks + intermediates) adds more per section.

---

## Phase 2 — build the cohort manifest

This replaces "edit `[paths]` and hope". The manifest is the contract: donor set, exact image
ids, and a fingerprint of every input file.

```bash
phenocycler manifest template --out data/qupath_manifest_spec.json
# edit the spec: expected_donors, and one entry per image
phenocycler manifest create --spec data/qupath_manifest_spec.json \
                            --out  data/qupath_manifest.json
```

Spec fields per image: `donor_id`, `image_id` (the exact QuPath `Image` value), `qptiff`,
`cell_geojson`, optional `nucleus_geojson`, `pixel_size_um_x/y`, `panel_id`, optional
`channel_map`. `defaults` carries the shared `measurement_csv` and a `segmentation_version`
string — record the QuPath project commit or model version there; it is part of the run's
identity.

`expected_donors` is checked against the images, so a donor you meant to include but forgot to
add an image for is an error, not a silently smaller cohort. Remember the exclusions above:
listing `6457` or `6579` will raise.

---

## Phase 3 — smoke test before committing to the queue

All of this runs on a login or dev node in minutes.

```bash
pytest tests/ -q                                     # ~15 s, fully synthetic, no data needed
python -m phenocycler.pipeline status                # no manifest yet -> tells you so, exits 0
python -m phenocycler.integration.pipeline --status
python -m phenocycler.integration.manifest --report  # donor <-> Xenium pairing table

python -c "from phenocycler.marker_registry import load_registry
r = load_registry(); print(r.registry_version, len(r.markers), 'markers')"
```

Once the manifest exists, `status` becomes the real thing: per-stage `CURRENT` / `STALE`, with
the reason a stage is stale. Run it before and after every job.

Then time one stage on real data before sizing the job:

```bash
/usr/bin/time -v python -m phenocycler.pipeline run --only ingest 2>&1 | tee logs/smoke.log
```

Read wall time and peak RSS out of that log. REDSEA's peak is driven by your image dimensions
(a full-resolution `int32` mask plus the distance transform's float64 + 2×int64 arrays over the
same `(H, W)`), so print one section's shape before believing any memory number.

---

## Phase 4 — run it

One job. Resource policy goes on the submission command, not in the script:

```bash
sbatch --account=<account> --qos=<qos> --cpus-per-task=8 --mem=120G \
       --time=2-00:00:00 scripts/slurm/run_full_pipeline.sh config.ini
```

The script `cd`s to the repo root, reads `PHENOCYCLER_JOBS` or `SLURM_CPUS_PER_TASK`, and
`exec`s `python -m phenocycler.pipeline run`. It does **not** activate conda — do that in your
submission wrapper, or use `conda run -n phenocycler-analysis`.

### The eight stages

| stage | writes | what it does |
|---|---|---|
| `ingest` | `cells/` | DuckDB: QuPath CSV → partitioned parquet, area floor applied |
| `geometry` | `geometry_qc/` | segmentation QC — area, solidity, raster/polygon agreement, duplicates |
| `redsea` | `redsea/` | lateral spillover compensation (`comp_mode=1`, `norm_form=donor`) |
| `expression` | `expression/` | which markers are measurable per donor |
| `controls` | `reference_controls/` | retained reference-control diagnostics (not RESTORE calls) |
| `calibrate` | `marker_evidence/` | RESTORE donor × target/reference fit → tri-state evidence |
| `type` | `assignments/` | RESTORE ordered residual → broad identity |
| `states` | `cell_states/` | process states, kept separate from identity |

Plus `qupath_export/` (UUID-keyed CSVs for round-trip), `audit/` (the models themselves) and
`manifests/` (immutable stage contracts).

Useful subsets:

```bash
python -m phenocycler.pipeline run --through calibrate   # dependency prefix
python -m phenocycler.pipeline run --only type states    # prerequisites must be current
python -m phenocycler.pipeline run --no-export           # skip the QuPath handoff
python -m phenocycler.pipeline export                    # export a current result later
```

Runs are **content-addressed**: `run_id` hashes the cohort manifest, the marker registry, the
ordered residual rules, the scientific configuration and the source tree. Change any of them and you get
a new run rather than a silently mixed one. That also means a code change — including this
merge — changes every `run_id`.

### Reading the result

`status` is the check, not the job's exit code. A stage reports `STALE` with the reason when
its inputs or method version moved.

The important honesty in `assignments/`: a cell is **Unavailable** when geometry
fails or a higher-priority gate cannot be evaluated. That is a real answer, not
a negative and not permission to try a lower gate. `Other` is reserved for an
explicitly all-negative marker profile.

---

## Phase 4b — what "positive" now means

RESTORE is donor-local. For each active target/reference pair, it fits the two
mutually exclusive groups with NNMF using raw QuPath `Membrane+Cell` features
where available (validated `Cell` fallback otherwise). The ordered
reference-positive/target-negative group defines the target background. The
production divisor is its `p99` under `control_definition=reference_only`.

A measurable analysis-eligible cell is `positive` exactly when its REDSEA
whole-cell target value is greater than that donor-marker divisor; it is
`negative` at or below the divisor. A missing measurement, geometry exclusion,
or technical pair failure is `unavailable`, never a negative fallback.

RESTORE also stores `restore_normalized = REDSEA value / divisor`, the active
reference, fit stability, control counts, feature set, and divisor in
`audit/restore_pair_models.parquet`. Re-running one donor never changes another
donor's fit.

### The knobs, and where they live

Marker behaviour is **data, not config**: `phenocycler/marker_registry.json` holds each
marker's RESTORE reference, family, and activation status, and
`phenocycler/ordered_residual_rules.json` holds the ordered compartment gates. Editing either changes `run_id`,
which is the intended way to make a change traceable.

`config.ini` `[restore]` records the frozen control definition, divisor
statistic, feature compartments, technical-failure policy, and exhaustive raw
compartment input. It also fingerprints the accepted donor 6566 INS
`CONFIRMED_ABSENT` review in `restore_expression_reviews.csv`. `[redsea]` and
`[geometry_qc]` hold their operator knobs.
`redsea_norm_form` is **validated**: anything other than `donor` raises, because the
mass-conserving form is what makes the compensation defensible.

### The taxonomy

The production type stage applies the historical hierarchy as an ordered
residual, not as five simultaneous competitors:

| priority | compartment | positive markers |
|---:|---|---|
| 1 | `Immune` | CD3e, CD68, CD11b, CD79a, CD11c, Iba1, CD206, MPO, CD4, CD8 |
| 2 | `Endothelial` | CD31, CD34, Podoplanin |
| 3 | `Epithelial` | E_cadherin, Pan_Cytokeratin, Ker8_18, EpCAM, Keratin_5, INS, GCG, SST |
| 4 | `Neural` | B3TUBB |
| 5 | `Mesenchymal` | Vimentin, SMA |
| 6 | `Other` | all 24 markers explicitly negative |

The first positive gate claims the remaining cell population. A gate with no
positive but any unavailable marker blocks every lower gate and emits
`Unavailable`. This default is broad-only; compatibility-specific columns
repeat the broad label. The probabilistic subtype implementation runs only
when an exact signed bundle is explicitly configured.

`docs/METHOD_HISTORY.md` records why the attempted migration was reversed, why
the ordered residual hierarchy was *not* restored with it, and which later
safety improvements were retained.

---

## Phase 5 — integrating with Xenium

Second environment, second orchestrator.

```bash
conda env create -f environment-integration.yml
conda activate phenocycler_integration && pip install -e '.[integration]'

python -m phenocycler.integration.manifest --report          # look before writing
python -m phenocycler.integration.manifest --build \
       --xenium-processed /path/to/Xenium_Analysis/data/processed
```

Build the manifest **yourself, first**. The `manifest` *stage* inside the integration pipeline
runs without `--xenium-processed`, so `xenium_h5ad`/`xenium_zarr` stay blank and
`import_xenium` falls back to the raw bundle — coordinates and morphology, **no cell-type
labels**. The run completes, the crosswalk is empty, and donor concordance compares nothing.

Two data-entry prerequisites the report will show you: `donor_overrides.csv` has two rows with
a blank `donor_id`, and `xenium_paths.csv` carries absolute `/orange/brusko/...` paths for 26
donors of which your cohort is 20 (set `[integration] xenium_root` to rewrite the prefix).

```bash
sbatch scripts/slurm/06_integration.sh                 # whole cohort, both tissues
ROI=pln_2 sbatch scripts/slurm/06_integration.sh       # one ROI
```

Read `data/integration/qc/integration_report.txt`, not the exit code: **PASS** = tissue Dice
≥ 0.80 and islet RMSE ≤ 150 µm; **WARN** = one gate; **FAIL** = neither. For lymph node, islet
RMSE is dropped from the decision rather than counted as a failure, since there are no islets
to measure against — and `structures`, `match` and `cellmatch` correctly report `n/a` there.

> ⚠️ **The PhenoCycler side of the integration is stale.** `export_pheno` still reads the
> retired `phenotype/broad/` and `restore_gated_redsea/` layout and has **not** been repointed
> at `assignments/` + `marker_evidence/`. The vocabulary crosswalk is likewise built against
> the old 7-class scheme. The config surface those 22 modules need is preserved so the package
> imports and the suite passes, but treat PhenoCycler→Xenium results as pending that rewire.
> The Xenium half is unaffected.

---

## Migrating an existing checkout

For the failed content-addressed run `c8412757d12de36ca056`, do not rebuild the
148 GB ingest. Use `pipeline adopt plan/apply` as shown in
`REFINEMENT_WORKFLOW.md`: the reviewed policy adopts ingest through controls
and the independent state stage, then recomputes RESTORE (`calibrate`), typing,
and QuPath export in a new run.

For a pre-manifest legacy checkout, build the cohort manifest first. Historical
RESTORE calls remain useful audit evidence but are not copied into a new run;
the new run refits RESTORE from fingerprinted raw compartments and REDSEA.

---

## Appendix — silent-failure checklist

```bash
# 1. Every donor in the manifest actually produced output
python -m phenocycler.pipeline status

# 2. Every donor x marker has an explicit RESTORE state and divisor disposition
ls "$PHENOCYCLER_DATA_DIR"/audit/restore_pair_models.parquet

# 3. Xenium arrived with labels, not coordinates only
grep -i "importing coordinates only\|donor_unknown" logs/integration_*.out

# 4. Registration verdicts, not the job exit code
cat "$PHENOCYCLER_DATA_DIR"/integration/qc/integration_report.txt

# 5. Crossmodal did not silently lose its anchors
grep -i "anchor" logs/integration_*.out
```

`crossmodal` refuses below `crossmodal_min_anchors = 8` verified protein↔gene anchors, and the
usable set is thinner than the 46 nominal 1:1 pairs because INS, GCG, SST and Vimentin are off
the Xenium panel. Lowering that number is a scientific concession, not a workaround.
