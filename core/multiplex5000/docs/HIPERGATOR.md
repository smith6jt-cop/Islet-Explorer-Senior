# Running the current upstream workflow on HiPerGator

This guide covers the presently supported PhenoCycler path on UF Research
Computing.

> **Current boundary:** no production target/reference pairs are verified.
> Ordinary execution stops after `expression`. Do not schedule `controls`,
> `calibrate`, `type`, or assignment export. The pair-independent `states`
> stage is explicit-only. Historical RESTORE,
> probabilistic, marker-debug, and Otsu commands are not cluster workflows.

## Supported stage graph

```text
QuPath measurement CSV + geometry + qptiff
  -> ingest
  -> geometry
  -> redsea
  -> expression
```

The cohort and stage manifests remain content addressed. Donor receipts allow
geometry, REDSEA, and expression work to overlap while preserving the same
cohort-level completion checks.

## Install

From a login node:

```bash
git clone <repository-url> Islet-Explorer-Senior
cd Islet-Explorer-Senior/core/multiplex5000
bash setup.sh
conda activate phenocycler-analysis
```

> The `phenocycler-analysis` conda environment is the HiPerGator recipe. On the
> development workstation the pipeline runs from the repo `.venv`
> (`docs/WORKFLOW_2026-09-02_MAP.md` § 7).

`setup.sh` installs the package and verifies the marker registry plus current
stage surface. It does not validate typing rules or RESTORE assets.

Keep environments and small source files on backed-up project storage. Put
large qptiffs, measurement exports, geometry, and run outputs on the project
filesystem appropriate for the allocation; do not copy them into Git.

## Declare the cohort

Create and edit the small specification:

```bash
phenocycler manifest template --out cohort-spec.json
```

Then build the fingerprinted cohort manifest:

```bash
phenocycler manifest create \
  --spec cohort-spec.json \
  --out data/qupath_manifest.json
```

Use absolute cluster paths in the specification when inputs are outside the
repository. The manifest, not a Slurm array range or filename convention, owns
the donor and image set.

Each manifest entry represents a **section**, not merely a donor. The stable
`donor_id` schema field is the partition key and therefore carries the section
key defined by `phenocycler.sections`:

| image | manifest `donor_id` | biological donor | ROI |
|---|---|---|---|
| `6539_Scan1.er.qptiff` | `6539` | 6539 | `panc` |
| `6539pLN_Scan1.er.qptiff` | `6539pln` | 6539 | `pln_1` |
| `6539pLN2_Scan1.er.qptiff` | `6539pln2` | 6539 | `pln_2` |

`expected_donors` consequently contains section keys and must match the images
exactly. `cohort.donors()` returns their unique biological donors for metadata
and Xenium joins. Repeating a bare donor across sections or using an unknown
section token is rejected during manifest creation; otherwise two coordinate
frames could enter one partition and be masked or normalized together.

Donor 6457 is centrally excluded and makes manifest resolution fail, including
through a section key of that donor (`6457pln2`, say). Do not bypass that
exclusion in a submission script. 6579's was lifted on 2026-09-17.

## Configure resources

The active scientific settings are cohort, geometry, marker-source, and REDSEA
contracts. Operational concurrency is controlled by:

```ini
[compute]
n_jobs = 4
geometry_workers = 2
redsea_workers = 1
downstream_workers = 4
use_gpu = false
gpu_device = 0
```

`n_jobs` is the total task ceiling. Geometry, REDSEA, and expression use their
group limits within it. If GPU REDSEA is enabled, keep one REDSEA task per GPU.

No `[restore]` or `[typing]` setting authorizes production work in the current
workflow.

## Submit

The provided launcher invokes the ordinary pipeline command:

```bash
sbatch \
  --account=<account> \
  --qos=<qos> \
  --cpus-per-task=8 \
  --mem=120G \
  --time=2-00:00:00 \
  scripts/slurm/run_full_pipeline.sh config.ini
```

Despite the historical filename, the launcher now follows the current stage
surface and stops after expression. The equivalent interactive command is:

```bash
python -m phenocycler.pipeline run \
  --config config.ini \
  --pipelined \
  --through expression
```

Do not append `--through calibrate`, `--only type`, `--no-export`, or a later
export command.

## Monitor and resume

Check the content-backed state rather than inferring success from a directory:

```bash
python -m phenocycler.pipeline status --config config.ini
```

Expected current stages are:

```text
ingest
geometry
redsea
expression
```

Adoption, resume, and recovery may reuse only validated upstream artifacts.
Historical controls, marker evidence, assignments, state calls, and QuPath
exports must not be linked into the current run.

If a job stops, rerun the same command. Current manifests and donor receipts are
validated before reuse; a non-empty partition without its valid receipt fails
closed.

## Development diagnostics

The current pairwise review is deliberately outside `data/runs/` and uses only
two upstream comparisons:

- raw CD3e versus raw Ker8_18;
- exact QuPath `Iba1_norm` versus raw Ker8_18.

CD79a is excluded as a pair candidate because its expression is too sparse in
normal tissue.

Refresh the narrow QuPath projection and run the rank-only overlap audit:

```bash
python scripts/review/qupath_fundamental_inputs.py \
  --donors 6436 6450 6566

python scripts/review/reference_overlap_audit.py \
  --donors 6436 6450 6566
```

These outputs are mutable development artifacts. They create no divisor,
positive call, cell type, stage manifest, or production claim. The mixed
Iba1/Ker8_18 comparison is evaluated by donor-local ranks, never by a shared
numeric threshold.

## Integration boundary

The Xenium integration environment and its data caveats remain documented in
[`INTEGRATION.md`](INTEGRATION.md). The current upstream PhenoCycler workflow
does not produce assignments for integration. Do not use a historical
`restore_gated_redsea` or phenotype directory as a substitute.
