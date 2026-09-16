# DeepCellTypes pilot — provenance and pins

Cell-type inference for the Phenocycler cohort using **DeepCellTypes**, a pretrained
panel-agnostic vision-language model (raw pixels + label mask + channel *names* → cell type).
It derives no per-marker threshold and needs no target/reference pair, so it does not revive the
retired RESTORE lane and does not depend on the blocked `controls`/`calibrate`/`type` stages.

Nothing here writes a `StageManifest` (no current stage writes one) or mutates the stable
`data/{cells,geometry_qc,redsea,expression}` roots or the historical `data/runs/` tree.
`pipeline.py` computes no `run_id` at all — the content-addressed run framework it fed was
removed on 2026-09-16 — so nothing under `scripts/**` can re-key a run that no longer exists.
Pilot output goes to the mutable repo-root `data/development/`.

## Software pins

Pinned **by commit**, declared in the `dct` optional-dependency group of
`core/phenocycler/pyproject.toml`:

| Package | Repo | Commit |
|---|---|---|
| `deepcell-types` | https://github.com/vanvalenlab/deepcell-types | `59d74a97711848dabc30fc01e242b6b1e4d3d3dc` |
| `deepcell-auth` | https://github.com/vanvalenlab/deepcell-auth | `523a4dd5d29941e076ffd378c6b032122589bdb2` |

`59d74a9` was upstream `master` on 2026-08-21. **The version string cannot substitute for the
commit:** `deepcell_types.__version__` is `'0.1.0'` on every master commit, so an unpinned
reinstall can change `resolve_supported_marker`, the crop/normalization path, or the abstention
default with nothing in the output able to show it.

`deepcell-types` requires Python `>=3.11`, above this package's own `>=3.10` floor. The repo
`.venv` (3.12.3, torch 2.13.0+cu130, CUDA on an RTX 4000 Ada / 20 GB) is the supported target.

Both commits are asserted against the installed `direct_url.json` files by
`tests/test_deepcell_crosswalk.py`, so a silent reinstall turns a test red instead of turning
results wrong.

## Model checkpoints — NOT installed by pip

`pip install -e 'core/phenocycler[dct]'` installs the code, not the weights. Weights come from
the authenticated `users.deepcell.org` portal through `deepcell_auth`, and land in
`~/.deepcell/models/`:

```bash
export DEEPCELL_ACCESS_TOKEN=...        # already set in ~/.bashrc
python -c "from deepcell_types.utils import download_model; print(download_model())"
```

`deepcell_auth/asset_manifest.yaml` is the single source of truth for keys and hashes. Read them
from `deepcell_auth.load_manifest()['models']['deepcell-types'][version]` at runtime rather than
transcribing them; reproduced here only so a copy obtained by any other route can be verified:

| Version | Asset key | md5 |
|---|---|---|
| `2026-06-15` (default; from-scratch resMLP, paper Fig-3c headline) | `models/deepcell-types_2026-06-15_resmlp.pt` | `b819a7e0b177ad5330394eab3c6c7ad8` |
| `2026-06-23-ptft` (SSL pretrain→finetune arm) | `models/deepcell-types_2026-06-23_resmlp_ptft.pt` | `402e94c103c5e489a57433cd107009d3` |

The 2025 CLIP checkpoints (`2025-06-09`, `2025-06-09_public-data-only`) predate the current
architecture; `download_model()` refuses them rather than handing back something that fails at
`load_state_dict`. Do not revive them.

**Pass an absolute `.pt` path as `model_name`, never a version string.** `_resolve_model_file`
(`predict.py:227`) routes a known version through `download_model()` on *every call*.

### Portal status

**Resolved 2026-08-23. Both checkpoints are in `~/.deepcell/models/`** and every result in this
lane was produced from them:

```text
deepcell-types_2026-06-15_resmlp.pt        20,849,567 B
deepcell-types_2026-06-23_resmlp_ptft.pt   20,850,427 B
```

*History, kept because it cost a day and would otherwise be re-diagnosed:* on 2026-08-21
`users.deepcell.org` returned nginx **503 on every path**, including `/api/getData/`, and served an
acme-companion self-signed fallback cert instead of its Let's Encrypt cert. That was verified
upstream, not local — `deepcell.org` and `github.com` presented valid CA-issued certs from the same
host, and `-k` did not help because the backend itself was 503. No local configuration change could
have worked around it, and none was needed: the portal recovered on its own.

## Inference needs no TissueNet archive

The model resolves channels against the packaged `deepcell_types/vocab.json` — **278**
standardized markers, **51** cell types — which arrives with the pip install. The ~2.4 TB
`expanded-tissuenet.zarr` release archive is required only for training and for two of the
paper's figure notebooks. It is not part of this pilot.

## The figure-reproduction repo is not a weights source

`literature/external/DeepCellTypes-2026_Wang_et_al` (upstream
https://github.com/vanvalenlab/DeepCellTypes-2026_Wang_et_al, commit `4a06e97`, retrieved
2026-08-20, Modified Apache License 2.0 with a non-commercial/academic carve-out) reproduces the
paper's figures. It has no torch dependency, no model code, and no checkpoint. Its
`data/required_datasets.yaml` lists 20 *derived outputs* (prediction CSVs, a CLS-embedding NPZ,
metrics JSONs) on the `deepcelltypes-2024-wang-et-al` S3 bucket, which denies anonymous access as
of 2026-08-21. None of it is needed to phenotype this cohort.

Two of its committed files *are* useful, and need no download:

- `notebooks/dct_figures/scoring.py` — the paper's own hierarchical-collapse scoring
  (`CELL_TYPE_HIERARCHY`, the ordered 51-class `CT2IDX`, `adjust_conf_mat_hierarchy`,
  `prf_from_cm`, `score_csv`), pure numpy/pandas, vendored out of the training code.
- `data/lineage_mapping.yaml` — cell type → lineage, the seed for `broad_class_crosswalk.json`.
  It covers only **44 of the 51** vocabulary classes, so it cannot be adopted verbatim.

Assert this checkout's HEAD matches `4a06e97` before emitting any `lineage` value.

The checkout lives at `literature/external/DeepCellTypes-2026_Wang_et_al`, which `.gitignore:24`
(`literature/*`) already covers. It was moved there from the repository root on 2026-08-24: as an
embedded git repository at the root it was neither tracked nor ignored, so `git add -A` would have
staged it as a **gitlink** — a submodule pointer, in a repo whose `CLAUDE.md` states there are no
submodules.

Its working tree is dirty — `notebooks/data_statistics.ipynb` and `notebooks/marker_positivity.ipynb`
were opened and run locally. That does not affect the HEAD assertion or the two files above, but it
means `git -C <checkout> status` is not clean and should not be treated as a corruption signal.

## Files here

| File | Role |
|---|---|
| `deepcell_types_crosswalk.json` | Phenocycler OME `channel_name` → DCT marker name. Holds **only** names the installed package cannot resolve itself, so it never silently forks from upstream. |
| `broad_class_crosswalk.json` | DCT class → broad compartment, **and the label space**: `in_label_space`, `exclusion_reason` and a written reason for every one of the 25 excluded classes. Read its `policy` block before changing anything. |
| `label_space.py` | Restricts the 51-way softmax to the 26 reachable classes. Exactly a renormalization over the kept columns — algebraically identical to masking the excluded logits. |
| `anchors.json` | The 28 gate anchors, **generated** by `concordance_gate.derive_anchors()` and frozen. Hand-editing it is what the tests catch. |
| `concordance_gate.py` | The pre-registered gate: does the model agree with this workflow's own REDSEA measurements, better than with anything else? Reports a permutation null and the 0.920 plane-identity ceiling beside every statistic. |
| `regions.py` / `regions.json` | The frozen 8-donor × 8-region sampling frame. Content-blind hash draw, reproducible from the donor id alone. |
| `window_mask.py` | Windowed rasterization for the seven donors with no archived mask, reusing `redsea._draw_geometry` unchanged. Verified bit-identical over 32 windows on four donors. |
| `arms.py` / `run_factorial.py` | The confound factorial: one variable per arm over concentric windows around a fixed anchor, scoring only the common 1024-px core. |
| `gates.py` | Fail-closed `Verdict` helpers. Absent, non-finite or under-powered evidence is a FAILURE, never a pass. |
| `tiff_window.py` | Windowed reads from a tiled or stripped TIFF without zarr. |
| `plane_identity.py` | Gate: does each marker *name* point at the right image *plane*? |
| `feasibility_tile.py` / `figure_feasibility.py` | The 2026-08-24 single-tile check and its figure. Superseded as a decision procedure by the gate above; kept because its four arms are what the results appendix records. |
| `tests/` | Run with `python -m pytest scripts/dct/tests -q` from `core/phenocycler/`. Not in the default `testpaths`. |

## Where the decision lives

`core/phenocycler/docs/DEEPCELL_TYPES_REVIEW_2026-08-24.md` — the pre-registered protocol, its
bar, and the results appendix. It supersedes the 2026-08-21 document, whose body is left unedited
on purpose.

## Naming trap

DCT's `Langerhans` class is the **Langerhans cell** — an epidermal dendritic cell, lineage
Myeloid — **not** the islet of Langerhans. In a pancreas project this is a mislabel waiting to
happen. Islet endocrine classes are `AlphaCell`, `BetaCell` and the generic `Endocrine`.

## 2026-09-05 — Track E.3 executed

`gates.py`, `tiff_window.py`, `window_mask.py` and `label_space.py` (with their tests) now live in
`scripts/common/`; this lane imports them from there. `broad_class_crosswalk.json` stays here because it
is this lane's decision record. The review document's "discard `scripts/dct/`" instruction is now
executable without breaking `scripts/restain/tiles.py` or `scripts/celltype/usability.py` (whose dead
`require_evidence` import was removed).
