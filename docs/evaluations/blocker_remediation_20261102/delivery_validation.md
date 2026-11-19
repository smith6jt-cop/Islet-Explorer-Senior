# Delivery validation and independent review

The final delivery checker is `scripts/evaluation/verify_remediation_delivery.py`. It refuses an unfinished acquisition ledger, incomplete cohort evidence, missing inspections, changed evidence references, approval claims, and incomplete project history. A final pass is recorded only in `outputs/blocker_remediation_20261102/delivery_verification.json`; focused checks and synthetic tests alone are not a completed-cohort result.

Run after the actual full-cohort reports, history and documentation are complete:

```bash
PYTHONPATH=core/multiplex5000 .venv/bin/python scripts/evaluation/verify_remediation_delivery.py --build-manifest
PYTHONPATH=core/multiplex5000 .venv/bin/python scripts/evaluation/verify_remediation_delivery.py
```

The first command also creates `artifact_manifest.json` and `provenance.json` in this directory. The second checks the saved manifest against the current files. Current source files require their actual current SHA256; they cannot pass by substituting historical Git content. The separate existing frozen-baseline checker retains its explicitly scoped historical-source fallback and independently checks that baseline evidence remains unchanged.

The manifest distinguishes freshly hashed source/evidence files from scientific files bound to previously completed full-byte audits and exact receipts. It also distinguishes current full image hashes from historical sampled hashes and binds the 24 recovered v20 files to their exact content-recovery audit. These categories do not claim a second read of the large scientific files. Missing, duplicated, substituted or differently classified bindings fail verification. Execution producers remain the receipts and captured runtime sources; a later delivery commit is not substituted for the execution producer.

The 25-section protein cohort also has 15 absent optional `pln_blocks:<donor>` inputs on sections with zero embedded pLN cells. These remain distinct `declared_optional_absence` bindings with null hashes and sizes. They are accepted only against the exact original receipt, the optional source declaration bound to its launch snapshot, the expected Stage 05 path, and the hashed full-section identity audit showing zero pLN cells. Unknown or required absences, nonzero pLN, changed files, dangling symlinks and conflicting receipt identities fail. The generic audit's original unavailable rows are preserved; absence is never counted as a successful byte hash.

The checker covers the exact 201-acquisition roster, all 19 dedicated pLN and 25 pancreas protein candidates, all 65 Xenium morphology sources with the three recorded damage exclusions, 32 reviewed MALDI scaling figures, 12 real same-slide pair dispositions, 60 exact CODEX checks, and 24 v20 content bindings. It keeps synthetic matching and batch tests separate from real accepted geometry. It checks archived observation evidence, project SQLite/export counts and hashes, MCP read/retry/export, disallowed approval, and unchanged global history.

Independent review found an arbitrary-selection bug in the MALDI Xenium source guard: duplicate experiment receipt names could select the first matching entry. The guard now requires exactly one matching experiment binding. Regression tests reject both identical duplicates and conflicting duplicates. A separate review of all 65 actual inventory bindings verifies uniqueness and checks all 62 validated source hashes against the existing receipt-bound source report; it does not rerun the image audit or alter its historical producer. Evidence is in `outputs/blocker_remediation_20261102/final_guard_review/source_binding_review.json`.

The recorder's immutable evidence archive, explicit source-index relocations, and historical full-versus-sampled hash verification were also reviewed. The archive checks bytes before recording, preserves original reference identifiers, and refuses corrupt existing archives or paths outside the project. No additional concrete defect was found in those paths. Focused delivery, recording and hash-algorithm tests pass 47 checks, including 14 optional-absence cases. The MALDI source and geometry tests pass 13 checks. Their reports are under `outputs/blocker_remediation_20261102/final_guard_review/`. These test counts do not grant scientific approval or release activation.
