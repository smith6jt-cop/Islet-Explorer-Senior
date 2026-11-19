# CODEX app and DuckDB fixes: recipient review

Received and read app commit `f673fb2` and security commit `b5fae5c` against `899a163`. The [receipt](receipt.json) pins full commit IDs, source/test evidence and the documentation snapshots. **No port, scientific relabelling, image build or deployment was performed.** The earlier [portability delivery](../blocker_remediation_20261102/README.md) remains unchanged and all eight verification checks pass.

The source worktree advanced after the requested commits. This review remains bound to those commits. The supplied display note and two app evidence files were uncommitted at capture; immutable copies preserve exactly what was reviewed. Final visual checks and a new dev image are outside this receipt's verification.

## Reproduced portability defect

`macrophage_collapse_comp()` mishandles a genuine raw Macrophage category alongside marker-profile macrophages when no split marker is selected. It merges core proportions, but then drops the new `peri_prop`, `peri_count`, `enrich_z` and `min_dist` destination columns. It also overwrites the residual peri count before reconstructing enrichment from the original members. Counts of 10 and 20 with expected counts 5 and 10 require z=3.872983; the internal result is 3.651484 before it is dropped.

The fix must compute from immutable original member columns and retain every merged destination family. Preserving columns alone is insufficient. Regression coverage should verify core/peri proportions, counts, reconstructed enrichment and minimum distance for mixed residual/profile inputs. The supplied current cohort has only twelve macrophage profiles, so this finding does **not** establish a current-cohort failure. A separate aggregate-only donor-pair error reproduces in the parent commit and is pre-existing.

## Independent checks and limits

Seven changed R suites pass 273 expectations across 51 blocks; three Python synthetic trajectory tests pass. The app review also passes 157 receipt consistency/helper-hash checks. All 22 app files are bound to `f673fb2` and identical at `b5fae5c`. The supplied whole-cohort partition reconciles 3,199,606 cells, and the 48 scenario roster is exactly 16 macrophage selections across core/peri/both. Eight actual receipt-declared input/helper files match their hashes. The reported full 588-assertion suite and real-data scenario calculations were not rerun here.

The security commit's raw diff and seven resulting files exactly match the reviewed `378843e` change. All 33 isolated patch/guard checks pass. Fresh patching reproduces the reviewed archive and changes only the two declared decoder source files. The patch receipt alone is not complete runtime attestation: helper/archive hashes must be verified and the newly built image must pass its installed-decoder checks. Historical image evidence belongs to its historical producer.

## Shared-rule adaptation needs

Keep parent totals separate from disjoint additive leaves. Map residual and unresolved states through the receiving tissue/panel schema; CODEX display filtering must not erase MultiPlex5000 ambiguity or unavailable evidence. Preserve original source labels/counts and record display exclusions and denominators separately. Donor-composition object memberships are not unique cohort cells.

The slash-label correction reads hidden Conflict memberships directly from raw regional columns, even when the pair lookup omits them. Conflict counts do not enter No subtype. Displayed bars and exported proportions of original cells intentionally use different denominators.

Live DPT/PAGA filtering changes computational inputs and may change a newly fitted ordering. Preserve the feature list and fitting provenance separately from frozen UMAP/Leiden results and scientific source identities. Layout changes received static review only.

Detailed evidence: [app review](../../../outputs/codex_app_fixes_receipt_20261103/app_review.json), [security review](../../../outputs/codex_app_fixes_receipt_20261103/security_review.json), [immutable documentation and arithmetic checks](../../../outputs/codex_app_fixes_receipt_20261103/root_evidence_review.json), and [unchanged prior delivery](../../../outputs/codex_app_fixes_receipt_20261103/prior_delivery_verification.json).

The first two attempts to notify the CODEX task failed because its internal MCP HTTP transport was unavailable. This receipt records the review; it does not falsely claim successful upstream delivery.
