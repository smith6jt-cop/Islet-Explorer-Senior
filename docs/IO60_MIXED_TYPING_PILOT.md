# IO60 mixed pancreas and pLN typing pilot

`io60_mixed_protein_pilot` runs Stage 06 for donor 6450 through the same
`multiplex5000.stages.typing` route used by `mplx5k`, the graph, and the
notebooks. It reads a completed and accepted IO60 Stage 00 and Stage 05
release, rebuilds the selected cohort-wide marker scale, executes the
pancreas hierarchy from the installed package, and evaluates the full 6551
lymph-node phenotype sheet on pLN cells from the same PhenoCycler scan.

The pilot is deliberately bounded. Pancreas refinement is fit on 6450 only;
it does not claim parity with the selected 25-donor representation or validate
any cell class biologically. pLN sheet outcomes remain `unique_candidate`,
`ambiguous`, `no_match`, or `unavailable`, never a served final class. No release
is activated and its Stage 06 gate requires human acceptance.

## Site inputs

Install the package from a commit or tag and set these absolute paths in the
site file's `sources` mapping (`MPLX5K_SITE_FILE` or
`~/.config/multiplex5000/site.json`):

| key | recorded source |
|---|---|
| `io60_upstream_release` | IO60 candidate containing accepted Stage 00 and 05 |
| `typing_geometry_root` | frozen `revised_cells`, `crossed`, and `qupath_fundamentals` tables |
| `typing_export_root` | v3 `structures` and `typing_islets` GeoJSON exports |
| `qupath_manifest` | v3 QuPath image, pixel size, and segmentation record |
| `donor_acquisition_runs` | donor-to-acquisition run table used by the selected hierarchy |
| `typing_registry_log` | selected `subpanel_log.json` criteria |
| `typing_shifts` | shifts record paired with Stage 05 alignment |

Every Stage 05 donor partition is declared and bound for the cohort scale,
although the hierarchy runs only for the pilot donor. The selected QuPath
structure and islet export files and each 6450 geometry table are declared
separately. The receipt records each input's path, size, and modification time;
ordinary input files are also hashed. The gate lists
inputs that still resolve into the old v3 worktree; those payloads must be
copied and verified under retained release inputs before promotion.

```bash
python -m pip install './core/multiplex5000[typing]'
mplx5k stage run 06_typing --project io60_mixed_protein_pilot
mplx5k gate show 06_typing --project io60_mixed_protein_pilot
```

The project contains one Stage 06 rule, so the graph and Stage 06 notebook
call the same module. A site may supply its own project file with the same
route and a different pilot donor or release ID. Do not point it at the active
release: this is a candidate review workflow.

## Evidence from the first local pilot

The 6450 scratch run retained 480,774 segmented cells with 480,774 unique
acquisition-scoped cell IDs. The pancreas hierarchy covered exactly 398,349
eligible non-pLN cells. The pLN route retained all 9,889 pLN cells, including
8,841 eligible cells, with the recorded `same-block` designation. The
installed-wheel run launched from `/tmp` produced byte-identical scale,
pancreas, pLN candidate, and mixed-cell Parquet tables to the checkout run.
The gate was not accepted. These numbers verify routing and replay, not the
correctness of the biological typing.
