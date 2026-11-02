# CODEX validation repair

The corrected content-bound replay passes all 60 exact numerical assertions across all 20 CODEX images. It reconciles 3,199,606 eligible cells to 3,525,977 measured objects and verifies the scientific inputs and baseline files stayed unchanged during the run. Ten focused tests cover content drift, proof-anchor drift, missing/duplicate acquisitions, incompatible runtime/thread settings, lossless reference parsing and mid-run input changes.

The repair is confined to reusable validation code under `scripts/evaluation`. The main checkout has no positivity package; the validator explicitly verifies and imports the pinned external `df3fe1d` snapshot. No scientific port was merged. The pinned numerical assertion functions, original computations, old fixture and historical failures remain unchanged.

The fixture uses the recovered producer journal's content bindings instead of comparing file dates with commit dates. It enforces the evaluated Python/package/thread contract before importing numerical libraries. It preserves the producer's default CSV parsing for computational inputs, while reading expected CSV outputs with round-trip float decoding. These are distinct table views; no numerical tolerance or scientific cutoff changes.

All 24 previously audited v20 release bindings were recovered byte-for-byte into `v20_recovered/` and independently reverified with an explicit relocation manifest. The matching levels/identity copies came from the CODEX repository's `authoritative_run`; the missing structure binding came from its hash-identical `binding_candidate.json`. The source run and release pointer were not changed. The v20 auditor now requires a fresh output directory, supports an explicit recovery map, and returns nonzero on failed checks.

Evidence: `exact_replay/replay.json`, `exact_replay/verified_provenance.json`, `fixture_tests.log`, `v20_recovered/recovery.json`, and `v20_final_verification/summary.json`. `completion_summary.json` reports totals; `artifact_manifest.json` binds the completed evidence.

This repairs the replay fixture and the 24 audited v20 content bindings. It does not claim that the entire v20 execution-input closure has been recovered, or confer biological acceptance or release activation. Original producer package versions were not historically locked; the checked runtime is explicitly the verified replay runtime.

All evidence paths above are relative to `outputs/blocker_remediation_20261102/codex/`.

The subsequently [received reference-context isolation handoff](codex_received_handoff.md)
is separate source evidence: 17 payload files verified, reviewed and not applied.
Its source/fix pins and adaptation requirements do not replace the running
candidate's producer or the pinned replay described here.

The later [component/Ki67 handoff](codex_received_component_handoff.md) is likewise
received and not applied. Local verification covers 17 payload files and 214
regular source files; its supplied scientific results were not rerun here.
Release/object failures remain in the supplied reports, and remote delivery was
not independently verified. Its pins remain separate from this replay.
