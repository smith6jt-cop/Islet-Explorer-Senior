# MultiPlex5000 execution board

Snapshot: 2026-10-21, `main` at `6a0561d`. This is a short control record for the
plan of record at `~/.claude/plans/evaluate-the-local-implementation-wondrous-seal.md`
and [branch dispositions](DISPOSITIONS_2026-10-20.md), not a replacement for either.
Update this board when a phase exits, a stage gate is accepted, or the active release changes.
Do not advance a phase because its PR merely merged.

| Work | Observed state | Exit still required |
| --- | --- | --- |
| A, preserve | The ten bundles under `/data/io60-archive/code_bundles_20261019/` pass `git bundle verify`; all 20 entries across its bundle and untracked-source `SHA256SUMS` pass (rechecked 2026-10-21). | Recheck the HPG Xenium and KINTSUGI preservation items before any related cleanup. |
| B, recover code | Complete per `DISPOSITIONS_2026-10-20.md`; PRs #25–#28 and #30 are merged. | Preserve the recovery manifests and producer receipts during later moves. |
| C, reconcile branches | Most stale branches were tagged and removed; the `hpg` remote is gone. | External repos and remaining worktrees wait for the later archive and data checks. |
| D1–D2 | Package rename, site settings, identity and relationships merged in PRs #29 and #31–#33. | Keep the temporary `phenocycler` alias until its callers move. |
| D3 | Graph, notebook route, gates and receipts merged in PR #34. Stage 00 coverage merged in #36; DuckDB read repair merged in #37. PR #38 adds acquisition panels, runs and tissue condition. | Rebuild Stage 00 from the final D3 code and obtain manual acceptance of its gate. |
| D4–D11 | Not built as unified stages. | Each stage uses the same functions on both routes and ends in reviewed figures and an accepted gate. |
| E–H | The complete replay, pilot, full rebuild, app switch, code retirement and one-release data cleanup have not happened. | Follow the plan's order. No active-release switch or deletion based on a Stage 00 candidate. |

The only built stage in the project file is `00_inventory`; `automatic_gates` is empty.
The current `io60-candidate-1/00_inventory` output has a receipt and a passing technical
gate but **no `accepted.json`**. It predates PR #38. The prior output was moved aside,
not deleted. The active release has not been switched to this candidate.

## Immediate sequence

1. Finish review of PR #38 as a Stage 00 change. Its green CI proves the code passes
   its tests; it does not approve the source interpretations or the figures. Keep
   D4 work out of the same PR.
2. Merge only after the review items below are resolved or explicitly recorded as
   unavailable. Rerun Stage 00 once from the merged commit into the candidate,
   preserving the old output and its receipt.
3. Compare the new `acquisitions.csv`, `coverage.csv`, `runs.csv`, both figures,
   source digests and `gate.json` with the recorded sources. The maintainer accepts
   the gate only after this review. Then start D4, the October 14 protein signal
   chain. Stages 09–10 and H1–H6 inference remain refused until processing is approved.

## Stage 00 review items

- **Missing is not zero.** At PR #38's tip, coverage formats a pLN designation
  with `if stated else ""`. The notes reader distinguishes `0` from `None`, but this
  expression renders both blank. Restore a visible zero for an explicit `no pLN`;
  leave an unstated block unknown. Test both cases.
- **IF acquisition versus export.** The registration inventory lists IF scans for
  6380 and 6457, while the maintainer's IF note says both had excessive punctate
  background and were not exported. Stage 00 currently counts each scan and says
  only that channel roles are unrecorded. Record `scan present; export skipped`
  with its reason in curated, tracked metadata and show that status in the gate.
  Do not invent channel roles or ingest the entire personal notes file.
- **Instrument-record provenance.** The 81 locally copied Xenium record files
  matched the sizes and SHA-256 values in `/data/xenium_records/MANIFEST.tsv` on
  2026-10-21. PR #38 reads and receipts the manifest and files, but its gate does
  not compare those claimed source digests to the copies. Make that a blocking
  check before accepting the Stage 00 gate.
- **Required fields.** PR #38 says each PhenoCycler slide ID is read. Its qptiff
  parser currently requires experiment, experiment date, scan time and markers,
  but accepts an empty slide ID. Require that field and test its absence.
- **Recorded relationships.** A section's block and a renamed Xenium region
  need a specific recorded source. A cohort-level "mostly different block"
  statement cannot fill an unknown section. The uncommitted `_blocks` edit in
  the PR #38 worktree currently assigns `different block` when a donor note is
  absent or the Xenium design is not sequential; that fallback must not land.
  Keep unresolved relationships visibly unknown and blocking where a later
  join would depend on them.
- **Figure review.** Check every pancreas's same-block and different-block
  pLN designation, the run partners and dates, panel labels, the six known
  Xenium folder/record name differences, the two tissue-condition rows, the
  five raw files absent locally, and each reported caution. A passing check
  does not turn a reported absence into a measured zero.

## Control rules

One integrator owns an implementation branch. Independent review may inspect
it but does not edit its worktree. A stage PR states its inputs, outputs, figures,
checks and any scientific decision it embeds. A receipt records what ran; only
an accepted gate allows the next stage. Keep `duckdb<2` until a separate CI run
passes with the bound removed. No active-pointer change, publication, archive
or deletion is implied by a passing Stage 00 gate.
