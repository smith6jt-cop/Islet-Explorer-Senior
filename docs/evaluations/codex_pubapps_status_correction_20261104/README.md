# PubApps status correction receipt

Commit `44daf5155ea14d015882aa876da1189be94dad5a` closes the four stale November 3 **status summaries** identified in our earlier review. The project README, `CLAUDE.md`, UFRC handoff and base-hardening guide now identify production release `pubapps-2026-11-ui-refresh-1659c9d` and the stopped local preview. The deployed-release receipt is byte-identical to the earlier package. [Previous receipt](../codex_pubapps_refresh_receipt_20261104/receipt.json) · [Correction receipt](receipt.json)

The workstation mirror and HPG publication have identical inventories and SHA-256 values for all ten files, including `SHA256SUMS`. Its digest is `a76c6fe35975cd50acd39afa8e719440185730110a639b92188360a1c02fb505`. The bundle passed `git fsck --full` with the acknowledged `a2f945a` base. Both email patches replayed directly in a disposable index and reproduced the upstream trees, ending at `eedc990814d9fab34ed5d6ca517129b650f8f330`. The two commits cover four Markdown documents and one handoff record. [Package verification](../../../outputs/codex_pubapps_status_correction_20261104/package_read.json) · [Git verification](../../../outputs/codex_pubapps_status_correction_20261104/git_review.json)

The review found three issues.

1. **Prior receipt hash.** `evidence/pubapps_refresh_agent_handoff_20261104.json` records a SHA-256 ending in `4434`. The receipt at recipient commit `1ce6007` and its current bytes both hash to `9404460d11edb025f60f6b4277706ac9e128f9c603374e4fdf3ede7f184b4438`. The difference is the final digit.
2. **Standalone README.** The package index lists four standalone corrected documents. Its own 900-byte `README.md` occupies the path for the 27,514-byte project README. The other three standalone documents match their tip blobs. The bundle and patch contain the corrected project README, SHA-256 `4099324befb941e684d2ef83e198accdb5fe04f1adc5ebe74fbacf081f1dbc63`.
3. **Lower README settings.** At the reviewed tip, the Development section still gives the requested loopback binding `127.0.0.1:3838` at lines 424–427 and labels schema `1` as current at line 456. The dated production receipt gives the deployed binding `10.13.41.21:3838` and schema `5`. [Status passage review](../../../outputs/codex_pubapps_status_correction_20261104/status_review.json)

The recipient disposition is **review only, patch unapplied**.
