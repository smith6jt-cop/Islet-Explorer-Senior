# Senior phenotyping steps 1–5 — pre-submodule originals (archived 2026-07-08)

These are the loose `scripts/senior/*.py` typing scripts that built the Senior 8-class broad-lineage
data **before** the pipeline was re-homed into the **`Phenocycler_Analysis` git submodule** — now the
single source of truth for steps 1–5 (cells → REDSEA → RESTORE (+extra) → hormone floor → broad
lineage → QuPath → figures). Kept for provenance only; nothing live imports or runs them.

Superseded by (submodule module):

| archived script | now |
|---|---|
| `build_cells_parquet.py` | `phenocycler.cells_parquet` |
| `redsea_full.py` | `phenocycler.redsea` |
| `restore_normalize.py` | `phenocycler.restore` (+ `--extra`) |
| `apply_hormone_floor.py` | `phenocycler.hormone_floor` (pipeline stage) |
| `assign_broad_lineage.py` | `phenocycler.lineage` |
| `marker_taxonomy.py` | `phenocycler.marker_taxonomy` |
| `plot_celltype_markers.py` | `phenocycler.figures` |
| `export_broad_class_for_qupath.py` | `phenocycler.qupath_export` |
| `redsea_reassess_diag.py` | `phenocycler.reassess_diag` |
| `redsea_full_validate.py`, `redsea_cohort_validate.py` | REDSEA validation harnesses; superseded by `reassess_diag` |

Parity was proven **byte-identical** (`assign_broad_lineage.py` vs `phenocycler.lineage`, donor 6539,
503,962 cells, 0 mismatches) before archiving. Run the pipeline now via:

```bash
pip install -e Phenocycler_Analysis          # once, inside `conda activate scvi-env`
bash scripts/senior/run_pipeline.sh --status
bash scripts/senior/run_pipeline.sh --only hormone_floor lineage figures --jobs 8
```

Downstream Senior scripts that used `from marker_taxonomy import EXCLUDED` were repointed to
`from phenocycler.marker_taxonomy import EXCLUDED`.
