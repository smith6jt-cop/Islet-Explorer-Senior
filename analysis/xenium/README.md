# `analysis/xenium` — Xenium 5K spatial transcriptomics

Absorbed from [`smith6jt-cop/Xenium_Analysis`](https://github.com/smith6jt-cop/Xenium_Analysis)
at commit `8062e554` (consolidation Phase 3). This is the **Xenium half** of the
PhenoCycler ↔ Xenium work; the protein half and the integration layer live in
`core/phenocycler/`, and the panel taxonomy lives in `core/panels/`.

The durable architectural guide — including the numbered do-not-regress decisions #1–#17 that
other documents cite by number — is [`CLAUDE.md`](CLAUDE.md).

## What runs here, and in what

**`~/miniconda3/envs/scvi-env`.** Verified 2026-08-28: scvi-tools 1.4.0.post1, scanpy 1.11.5,
anndata 0.12.4, torch 2.9.0+cu128 with CUDA available on the RTX 4000 Ada, and every API the
typer calls (`SCANVI.from_scvi_model`, `SCVI.train(load_sparse_tensor=)`,
`SCANVI.train(n_samples_per_label=)`, `sc.tl.score_genes(layer=)`).

**Not `xenium_analysis`.** That env is broken on this machine — `import scanpy` and
`import scvi` both raise `ImportError: cannot import name '_spropack' from
'scipy.sparse.linalg._propack'` from a mixed scipy/numpy install. Its `environment.yml` was
dropped on absorption rather than carried across as a recipe for a broken env.
`docs/installed_freeze.yml` is kept as the provenance record of what HiPerGator had.

**`rapids_singlecell` is absent and is not needed.** `try_import_rsc()` returns `None` with a
warning and **every call site has a complete CPU branch**. Decision #7's `LD_PRELOAD` recipe
is HiPerGator-only and now lives in `hpc/`.

## Two roots, deliberately separate

```python
ROOT      = Path(__file__).resolve().parents[N]   # analysis/xenium — CODE + tracked tables
REPO      = ...                                    # monorepo root
DATA_ROOT = os.environ.get("XENIUM_DATA_ROOT", "/data/xenium")   # the h5ads
```

Collapsing them — as the single HiPerGator `ROOT` did — would put `data/processed/` under
`analysis/xenium/`: 60+ GB of h5ads inside the repo tree, in a directory the repo-root
`.gitignore` pattern `/data/` **does not** cover, because that pattern is anchored.

## The panel taxonomy is generated, not copied

`core/panels` **is** XeniumPanelExplorer, imported as a subtree in consolidation Phase 2.
Everything `scripts/fetch_subpanels.py` used to fetch over HTTP from the pinned
`XeniumPanelExplorer@076d1192` is already tracked there — verified byte-identical on
2026-08-28 for `panel_roles.csv`, `panel_gene_roles.csv` (1,613,132 bytes) and the identity
cores. So the script reads locally:

```bash
python analysis/xenium/scripts/fetch_subpanels.py      # -> 129 files, 4 phases, 0 failed
```

One authority (`core/panels`, tracked, in git, pinned by decision #17); one derived cache
(`data/development/xenium_panels/subpanels/`, regenerable, gitignored, never hand-edited),
stamped with `SOURCE.txt` recording `git rev-parse HEAD:core/panels/data` — the tree object
id of the taxonomy itself, which moves if and only if a byte of it moves. `tests/` compares
the stamp against the live tree, so a stale snapshot fails a test.

## The typer

```bash
conda run -n scvi-env python analysis/xenium/scripts/finalize_panel.py \
    --tissue pancreas --audit-csv analysis/xenium/panel_audit.csv \
    --subpanel-dir data/development/xenium_panels/subpanels \
    --paths /data/xenium/processed/*/panc/*_phenotyped.h5ad

conda run -n scvi-env python analysis/xenium/scripts/hierarchical_phenotype.py \
    --tissue pancreas --output slim --dl-workers 0 \
    --subpanel-dir data/development/xenium_panels/subpanels \
    --paths /data/xenium/processed/<serial>/panc/<serial>_phenotyped.h5ad
```

> **2026-09-02:** `hierarchical_phenotype.py`'s `main()` does **not** yet define `--paths`
> (only `--samples`, which resolves the flat `processed/{s}/{s}_phenotyped.h5ad`), and a missing
> input is skipped silently with exit 0. Until the fix lands
> (`docs/PLAN_2026-09-02_AUXILIARY_PACKAGES.md`, Track B.1) the second command above cannot reach
> the nested `<serial>/<roi>/` files the builder writes.

`finalize_panel.py` must run **first** (it stamps `identity_panel`). `--dl-workers 0` is
mandatory — any value > 0 deadlocks the scANVI finetune at epoch 1 with the GPU idle
(decision #16). `--output inplace` overwrites with **no backup**.

## What was dropped on absorption, and why

| dropped | reason |
|---|---|
| `config.ini` | Its `[phenotyping]` block mirrored **retired** panels — `Epsilon = GHRL` is the single-gene argmax that mislabelled ~150K cells and that decision #12 exists to prevent; `Acinar = AMY1A` is the 1-gene panel `pancreas.py` documents fixing. Root `CLAUDE.md`: there is exactly one `config.ini`, and it is `core/phenocycler/config.ini`. A reference mirror that contradicts the canonical method is worse than none. |
| `notebooks/06_xenium_phenocycler_integration.ipynb` | An **unrun skeleton**: placeholder `SAMPLE_NAME`, a segmentation mask stamped from discs at cell centroids, one global `phase_cross_correlation` shift. Superseded by `core/phenocycler/phenocycler/integration/register.py`, which has two registration tracks, a px→µm scale fitted from the data, and vertex-wise polygon warping. |
| `xenium_paths.csv` | Superseded by `core/phenocycler/data/integration/xenium_paths.csv`, the de-brusko'd site-neutral relative-path version. Importing this one would reintroduce `/orange/brusko/...`, which was removed to zero references. |
| `donor_metadata_panc.xlsx` | It is data. Now at `core/phenocycler/data/raw/`, untracked, where `config.py` already points. |
| `environment.yml` | Names the broken `xenium_analysis` env (above). |
| subtree `.github/workflows/`, `.pre-commit-config.yaml`, `setup.sh`, `requirements-ci.txt` | GitHub reads only the root `.github/workflows/`; precedent is commit `658b9d7`. |
| `scripts/program_panels.py` | Superseded by `state_panels.py` per decision #17. |
| `scripts/R/xenium_analysis.R` | Self-described demo; hardcodes the broken conda env; duplicates the Python path. |
| `data/reference/embedding_genes_*` | The retired `embedding_panel`, listed under "Retired" in decision #16. |

The **thymus** lane (`tissue_markers/thymus.py`, `panel_audit_thymus.csv`,
`ingest_thymus.py`, `thymic_zonation_analysis.py`, `docs/THYMUS_PIPELINE.md`) is kept but is
**out of scope** for this cohort: `finalize_panel.py` and `_panel_roles.py` are
tissue-parametric, `core/panels/data/tissues/thymus/` is already in the subtree, and removing
it would mean editing shared code for no gain.

## `hpc/` is HiPerGator-only

The SLURM lane, `run_pipeline.sh`, `env_rapids.sh` and the benchmark env builders. They
parameterise on `XENIUM_ANALYSIS_ROOT` / `XENIUM_DATA_ROOT` and keep the HiPerGator paths only
as documented `${VAR:-default}` fallbacks. **None of it runs here** — `ssh hpg.rc.ufl.edu`
currently returns `Permission denied (publickey,keyboard-interactive)`.

## `scripts/benchmark/` is dormant

BANKSY · UTAG · KODAMA · MISTy on one identical input, with `eval_lib/` and a Tier-C
scorecard. It is scoped to spatial **domains**, on samples this cohort does not contain, and
neither of its two environments (`bench_py` conda, HiPerGator `R_libs`) exists here.
`eval_lib.{agreement,bio_recovery,domain_metrics,figures}` do import and run in the repo
`.venv` unchanged — the metric library is path-free, and it is the reusable half.

## Tests

```bash
cd analysis/xenium && ../../.venv/bin/python -m pytest -q tests/     # 77 passed
```
