# `scripts/xenium_typing/` — the Xenium producer comparison (request 2)

Development only: writes `data/development/xenium_typing/`, never the bases in `/data/xenium/processed`
(the incumbent typer's sidecars are read, not rewritten). Rules: root `CLAUDE.md`; `analysis/xenium/CLAUDE.md`
decisions #1–#17 (identity → mask, subclass → fine, state → score; no HVG selection; no library-size
covariate; `random_state 0`; `dl_num_workers 0`). Protocol: the plan's Phase 2 (`~/.claude/plans/make-a-plan-to-goofy-pnueli.md`).

| module | env | what |
|---|---|---|
| `pins.py` | `.venv` | donors, serial/run maps, base/sidecar/bundle paths, identity cores, panel snapshot stamp, `OUT` |
| `gate.py` | `.venv` | the one record shape (`cell_id, lineage_native, celltype_fine, confidence, reason, cycling` + `state_*_score`), `write_calls`, the label-free gate (`judge`): null-gated lineage award on the labels, spatial coherence vs within-section permutation, stability across replicates, cross-arm agreement, the same-section IF INS/GCG rank association (descriptive) |
| `nullgate_typer.py` | `.venv` | **X1** `scanvi-per-donor`: the incumbent sidecars read into the record shape (stochastic, no replicates → stability unjudged) |
| `joint_model.py` | `scvi-env` | **X2** `joint-run-aware`: one scVI over a 50k-cell hash subsample per donor, `batch_key = xenium_run`, Leiden × 5 seeds at 3 resolutions, null-gated cluster annotation, resolution by seed stability, per-cell kNN vote in the joint latent, fine tier by identity-core `subtype_hint` panels, cycling rule, state scores after typing |
| `baseline_analyzer.py` | `.venv` | **X3** `xenium-analyzer-baseline`: the instrument's `graphclust` (primary) and k-means ladder (replicates), annotated by the same null gate; `Unavailable(not_transferred)` until the Globus copy lands |
| `release.py` | `.venv` | the B.2 release for `phenospace.load_release`: every bundle cell, QC-dropped = `Unavailable(qc_dropped:<set>)`, non-null `evidence_positive`, `feat__if_INS/GCG` |

Run order: arms → `gate.judge(arm, donor)` (writes `gate/<arm>/donor_id=<d>.json`) → `release.py --arm <producer>`.
Producer choice follows the protocol's mechanical rule (verified lineages on both donors; then margin; then
cohort coverage).

`pins.donors()` reads the acquired cohort from `core/phenocycler/data/qupath_manifest.json`,
not the older build specification, and normalizes `0112` to `112` before checking central
eligibility. The 22 donors span 12 Xenium runs. Bundle selection stays keyed to
`xenium_paths.csv`: 6533 uses `0041393` and 112 uses `0041396`, both in MM036. Rows marked
`PhenoXeno` describe separate same-slide acquisitions and are excluded from this lane's
maps. Duplicate cohort rows still raise, including in `run_map`.

Tests include live cohort readback plus synthetic typing, missing-input, exclusion and
bundle-selection checks. Run from the repository root with the shared data mounted:

```bash
.venv/bin/python -m pytest -q core/phenocycler/scripts/dct/tests core/phenocycler/scripts/xenium_typing/tests
```
