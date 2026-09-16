# `scripts/report/` — the two-donor figure set (6436 T1D, 6442 ND)

Development lane (plan 2026-09-05, Phase 4). Static PNG (dpi 200) + PDF + a JSON fact sidecar per
figure under `data/development/report/<donor>/{pheno,xen,multi}/` and `_pair/` for the four
both-donor figures. `figure_manifest.json` lists every input with its sha256 (sampled for files
above 8 GiB, stated on the entry); `--verify` re-hashes and exits non-zero on drift.

```bash
cd core/phenocycler
../../.venv/bin/python scripts/report/run_report.py --preflight
../../.venv/bin/python scripts/report/run_report.py --section all            # or: pheno xen identity multi
../../.venv/bin/python scripts/report/run_report.py --section identity --only P09 X07 --donors 6436 6442
../../.venv/bin/python scripts/report/run_report.py --verify
../../.venv/bin/python -m pytest -q scripts/report/tests
```

## Figures

| section | ids | reads |
|---|---|---|
| `pheno` (multiplex proteomics QC) | P01 counts/eligibility, P02 areas + solidity, P03 eligibility reasons, P04 per-marker asinh densities with the Step-2 background mode, P05 REDSEA window (CD3e / Ker8_18) | the immutable run `30f6689e…` and UniFORM Step 2 only |
| `identity` PhenoCycler | P06 marker × label heatmap beside a label-permutation null, P07 spatial map, P08 diagnostic UMAP, P09 proportions by donor (pair; raw + evidence-positive), P10 proportions by QuPath `cell_region` (audit) and by identity-derived zone, P11 `ki67_donor_rank` per label (no cut) | `data/development/broadtype/release` |
| `xen` (spatial transcriptomics QC) | X01 per-cell QC, X02 drop reasons (overlapping sets), X03 negative controls + `metrics_summary` decoding metrics (pair) | bundle + base only |
| `identity` Xenium | X04 identity-core dotplot beside the nullgate award, X05 diagnostic UMAP (joint latent), X06 spatial map (QC-dropped = Unavailable), X07 proportions by donor (pair; raw, evidence-positive, low-signal share), X08 proportions by identity-derived islet zone (S3 structures; no QuPath annotation on Xenium), X09 cycling per label with the rule printed | `data/development/xenium_typing/release` |
| `multi` (sequential-slide) | M01 registration overlay + the serial-section verdict and the SAME-SECTION IF↔Xenium registration in separate boxes, M02 side-by-side common-lineage maps each in its own frame, M03 composition concordance over `comparable_for('pancreas')` with a within-donor permutation null (pair), M04 islet level (matched islets only where the verdict permits; else registration-free counts/sizes), M05 same-section IF vs RNA (per-gene Spearman ρ, IF hormone per Xenium label beside a permutation null, Xenium-labelled IF nucleus overlay), M06 niche-view concordance | releases + `data/integration/` stage outputs + `data/development/restain/` |

## Rules enforced here

- Nothing re-types a cell. Labels are read from the releases; `evidence_positive` is read as written.
- A missing release or stage output renders a **stated "not available"** panel with the reason; nothing is faked.
- Every proportion panel prints its denominator and is paired with its evidence-positive version.
- `Ambiguous` is hatched, `Unavailable` grey, and the common `Unassigned` class is never drawn.
- Axes start at 0 (asinh or symlog when the range demands it; never a percentile black point).
- Every statistic sits beside its null (label permutation, within-donor permutation, lineage permutation).
- Identity tier only: state markers appear in P06 as reported rows below the identity block, never as inputs; Xenium state subpanels are never drawn as identity evidence.
- The two registrations are never conflated (M01 boxes); the transform is drawn only in M01 while the serial-section verdict is FAIL.
- Title line 2 on every figure: development evidence, latch `NO_VERIFIED_TYPING` down, this method at this configuration on these donors. The footer prints run/content ids, code commit, and the release ids once they exist.
- Nothing is copied out of `data/development/restain/explorer_bundle/`. This was the DO_NOT_SEND
  hold; the maintainer lifted it on 2026-09-15 and the `FORBIDDEN_SOURCE` pin was retired with
  it, so this is now a statement of what the lane does rather than a prohibition.
