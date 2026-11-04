# November 4 PubApps incremental handoff

The HPG publication and workstation mirror match byte for byte. Both contain the seven listed payloads and `SHA256SUMS`, whose SHA-256 is `eaa9497a6b25c80a6fce24bec9c5a44ff8b5915c80e8f4e471e0170b4b2b27ef`. The named HPG directory resolves to `/blue/maigan/smith6jt/CODEX-handoff-20261104-pubapps-refresh.content`. The recipient disposition is **reviewed, adaptation unapplied**. [Machine-readable receipt](receipt.json) · [Package inventory](../../../outputs/codex_pubapps_refresh_receipt_20261104/package_read.json)

| Role | Git identity |
|---|---|
| Previously acknowledged handoff base | `5b72a9b022766f5ddb828845e5789bfcf3f15063` |
| Application, image and data producer | `1659c9dd01d8d409030cfb4618eb9f3edf56488e` |
| Parity validator | `7ffa7a1f38a95a23f311e62892df857f2f16b466` |
| Package tip with production records | `a2f945ac4c20d06a1354a38cfb93d5491f1f8cef` |

Two findings affect how to use the package.

1. **Use the bundle to recover the original history.** In a disposable repository containing the prerequisite base, the bundle passed `git fsck --full` and reproduced the upstream graph of 58 commits, including two merges. The email patch contains all 56 non-merge commits and matches a freshly generated upstream patch byte for byte. Direct linear replay stopped at message 27 (`ffeb729`) in `test-macrophage-labels.R`. Replaying the messages with `git apply --cached --3way --binary` from the base reproduced the exact tip tree `dcfe1f37028a00afb0d1939071ed9bda9c13ecc6`. [Git verification](../../../outputs/codex_pubapps_refresh_receipt_20261104/git_review.json)
2. **Four documents at the package tip still present November 3 status as current.** The affected passages are `README.md:27`, `CLAUDE.md:3`, `PUBAPPS_UFRC_HANDOFF.md:15` and `deploy/BASE_HARDENING.md:3`. For November 4 operations, use the dated production receipt, `docs/pubapps_ui_refresh_20261103.md`, `docs/workflow_user_guide.md` and `deploy/README.md`. Those documents identify the activated production release and record that the local preview was stopped. The workflow and deployment guides are available through the bundled Git history. [Guide review](../../../outputs/codex_pubapps_refresh_receipt_20261104/guide_review.json)

The three packaged evidence files match the Git blobs at the package tip. The production receipt records activation of `pubapps-2026-11-ui-refresh-1659c9d` at **2026-11-04 13:18:01 UTC**, with the previous release retained for rollback. The validator changed only `scripts/verify_release_parity.R`. The saved report has 25 distinct passes at a numeric tolerance of `1e-12`. The parity report's `release_manifest_sha256` identifies the data-root manifest (`fd21101e…`). The outer release manifest is `a0977425…`. [Evidence review](../../../outputs/codex_pubapps_refresh_receipt_20261104/evidence_review.json)

The access report has 15 distinct passes, comprising nine anonymous Basic Auth challenges and six authenticated requests, including a 16-byte TIFF range and WebSocket upgrade. Service health, host hashes and read-only mounts are producer observations. Our independent checks covered the package bytes, Git reconstruction and saved reports. The designated off-campus check, VM reboot test and confidential delivery of reviewer credentials are recorded as pending. The earlier 49-check browser run retains its original release identity.
