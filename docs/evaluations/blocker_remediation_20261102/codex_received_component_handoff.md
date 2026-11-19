# Received CODEX component and Ki67 handoff

The separate `CODEX-handoff-20261102-component-build` package was received,
verified locally and read. **It has not been applied.** The
[machine-readable receipt](codex_received_component_handoff.json) and
[file verification](../../../outputs/blocker_remediation_20261102/codex/component_handoff_local_verification_v1.json)
preserve its identity separately from the earlier
[isolation handoff](codex_received_handoff.md), the original CODEX parity baseline
and the active 25-section PhenoCycler run.

The independently computed SHA256 of `SHA256SUMS` matches the supplied value:
`7e091cea244305ceef160c4a2d83a7b519c20604f0ca46ede40481c84223b388`.
All **17 payload files, 2,056,899 bytes**, match that manifest. All **214 regular
source files** match the source manifest inside the archive. There are **238 tar
entries** when its 24 directories are included. No archive was extracted into a
recipient source tree. The bundle head and prerequisite header match source
`2150edd93ea379acd60cc36db3eb0b18eeb207d4` and history base
`29f45b35dc1566d4fce89d0e8813fa5b55743c4c`; Senior did not import or verify the
prerequisite Git history. Remote delivery is reported by the sender and was not
independently checked here.

The supplied scientific build is separately identified as
`b15f0fcb4566dcfaf8f7a02fcaa88c550e3f3beb`. Its reports claim 20 donors,
3,199,606 eligible cells, 110,721 preserved named v19 endocrine calls and 351
Lymphatic-only label changes in the mixed-control comparison. These differ from
the earlier initial-v20 diagnostic's 113,134 endocrine calls and 352 changes.
This receipt verifies package integrity and arithmetic within the supplied
summaries; it does **not** reproduce those scientific results.

The supplied reports retain **8,629 unresolved object-class labels, 2,472
component conflicts, 953 populated objects lacking a positive eligible member,
and three Single membership anomalies**. Their release preflight and strict
object checks fail. Assembly and export checks do not provide scientific
acceptance or resolve those failures.

The guidance read for a possible later adaptation is:

- Bind each selected component, calibration context and exact UUID membership
  explicitly; retain hormone-evidence replay and unresolved boundary records.
  Regenerate dependent outputs from the coherent cell assignments.
- Preserve binary `Ki67_pos` through cell, donor, object and distance exports,
  with per-image measurement/threshold provenance and CSV/Parquet parity checks.
  The sender reports 78,680 Ki67-positive cells. Ki67 remains a state and is not
  added to identity or the 31-marker embedding input.
- Keep fixed calibration references separate from live structure membership.
  A recipient port needs its own assay, tissue, schema and registry mapping;
  CODEX thresholds, selected v19 labels and cohort assumptions are not defaults
  for PhenoCycler, RNA, pLN or spleen.

No core code, current producer, scientific output, registry approval or release
activation changed. The earlier absence receipt remains a historical observation;
this new receipt records the subsequently delivered package without rewriting it.
