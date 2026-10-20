# Recovered 2026-10-08 producers — a record, not a runnable lane

These files were committed nowhere before 2026-10-20. They lived inside the v3 worktree's
`results/` tree (and one in a data folder), and they produced data that is current today:

- the protein measurements and calibration under
  `/data/io60-calibration/own_compartment_frozen_reference_repair_20261008` (the October 14
  chain, the canonical protein signal path), and the full-cell release and admission tables built
  on them;
- the conservative protein candidate and the reviewed pancreas and pLN RNA candidates the 10-14
  handoff selects;
- the H1–H6 evaluators and the RNA reader the cohesive lane snapshots and loads.

Each file is kept **byte for byte**. `RECOVERY.json` lists, for every one, its original path, its
sha256, and the receipts that name its path or record those exact bytes (28 of the 53 files in
this salvage have their sha256 recorded in a current receipt; the rest are named by path, by the
lane, or exist nowhere else).

**Do not run them from here.** They hard-code paths into the v3 worktree and `/data`, and some read
inputs whose own producers are still only in the v3 `results/` tree. They sit outside the live
producer trees (`core/phenocycler/scripts`, `scripts/senior`) on purpose: the producer guard in
`core/phenocycler/tests/current/test_no_superseded_work_is_reachable.py` polices live readers, and
these are not live. Phase D of the MultiPlex5000 plan ports what survives into
`core/multiplex5000` with declared inputs, and the replay checks there compare against the
generations these files made.

Everything else from the v3 worktree is preserved in
`/data/io60-archive/code_bundles_20261019/untracked_source/` (see its `MANIFEST.tsv`).
