# Fundamental-pair review — what the 2026-08-14 figures show

**Status:** development-mutable review evidence. `creates_cell_calls: false`,
`creates_thresholds: false`, `production_outputs_written: false`. This is rank/overlap
troubleshooting on upstream QuPath inputs. It creates no marker call, no threshold, no cell
type and no production pair, and it writes no run or stage manifest.

Written 2026-08-20. The figures were produced 2026-08-14 13:54 into
`data/development/fundamental_pair_figures/` and, until now, appeared in no document — the
`HANDOFF_2026-08-12_RESTORE_INTENSITY.md` "deliberately not completed" list still reads
*"Require a materially improved three-donor result before any 20-donor run"*, and the
three-donor result it asked for was sitting unread on disk.

## What was compared

Roster `fundamental-qupath-input-pairs-v4`, schema v3, method
`fundamental-pair-review-figures-v8`. Two comparisons, three donors (**6436, 6450, 6566**):

| target | measurement | reference | measurement |
|---|---|---|---|
| `CD3e` (T_cell) | raw QuPath `Cell: CD3e: Mean` | `Ker8_18` (Epithelial) | raw QuPath `Cell: Ker8-18: Mean` |
| `Iba1` (Myeloid) | exact signed `Cell: Iba1: Mean (norm)` | `Ker8_18` (Epithelial) | raw QuPath `Cell: Ker8-18: Mean` |

Mixed scales are deliberate and are evaluated **by rank only**; no shared numeric cutoff is
applied anywhere. CD79a is excluded as a pair candidate — too sparse in normal tissue.

Note on the Iba1 scale, because it is easy to misread: QuPath's norm is
`(raw − tissue p50) / (tissue p99 − tissue p50)`, so **1 is the tissue p99 anchor, not a
positivity threshold**, and legitimately-negative values are expected.

## Result

**`Iba1_norm` vs `Ker8_18` behaves like a fundamental pair.** Both correlations are negative
in all three donors and the joint upper tails are essentially disjoint:

| donor | Pearson r | Spearman ρ | q99 both |
|---|---|---|---|
| 6436 | −0.219 | −0.259 | 1 cell (0.01%) |
| 6450 | −0.179 | −0.188 | 6 cells (0.14%) |
| 6566 | −0.261 | −0.217 | 1 cell (0.00%) |

**`CD3e` vs `Ker8_18` does not.** Pearson r is ≈ 0, but Spearman ρ is **positive in all three
donors** — the two markers rank *together*, which is the opposite of what mutual exclusivity
requires of a target/reference pair:

| donor | Pearson r | Spearman ρ | q99 both |
|---|---|---|---|
| 6436 | −0.001 | **+0.114** | 99 cells (0.51%) |
| 6450 | −0.021 | **+0.146** | 4 cells (0.09%) |
| 6566 | −0.065 | **+0.268** | 10 cells (0.05%) |

The joint density shows the mechanism: CD3e sits at a mode of roughly 150–200 native units
even where `Ker8_18` is near zero, so it has **no resolvable negative population**. `Iba1_norm`
by contrast has a tight negative mode at ~0.0 with a sparse positive tail — which is what makes
it separable at all.

`rank_tail_overlap` agrees. Iba1/Ker8_18 is depleted relative to independent tails at every
quantile with no anomalies. CD3e/Ker8_18's depletion is shallower and donor-inconsistent
(enrichment ratio ~0.6 for 6436 against ~0.2 for 6450 and 6566 at the 10% tail), and its one
apparent enrichment spike — donor 6450 at the 0.1% tail — is **n = 1** and carries no
information. Below roughly the 1% tail all three donors fall to n = 0–2, so the strict tails
are uninformative for both comparisons; the usable range is 10%–1%.

## What this does and does not license

- It **does** say the two rostered comparisons are not equivalent, which the roster currently
  treats as though they were. On this evidence `Iba1_norm`/`Ker8_18` is the better-behaved of
  the two by a clear margin, and CD3e's positive rank correlation with an epithelial reference
  is a red flag that needs explaining before CD3e is used as a target against it.
- It **does not** verify a production pair. Three donors, rank-only, no thresholds, no calls,
  mutable workspace. `production_status` remains `development_no_verified_pairs`.
- Per the standing rules, none of this may be resolved by tuning toward an expected cell-type
  frequency, cross-donor uniformity or disease ordering. The per-donor spread — ρ ranging
  +0.114 to +0.268 for CD3e — is itself part of the finding, not noise to average away.

## Figures

`data/development/fundamental_pair_figures/` (mutable; may be refreshed in place):

- `fundamental_pair_joint_distributions.{png,pdf}` — 160-bin/axis joint densities on
  donor-native linear QuPath axes, with the r/ρ/q99 annotations quoted above.
- `rank_tail_overlap.{png,pdf}` — exact-rank upper-tail overlap and enrichment relative to
  independent tails, at the 0.9 / 0.95 / 0.99 / 0.995 / 0.999 quantiles.
- `crop_montage_6450.{png,pdf}` — raw-channel crops, deterministic rank-extreme selection,
  explicitly nonrepresentative by design.
- `figure_manifest.json`, `roster_used.json` — provenance and roster version.
