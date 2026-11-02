# Handoff: consolidation state, 2026-11-02

This snapshot updates the state recorded in the
[October 28 coordination handoff](HANDOFF_2026-10-28_SENIOR_CODEX_COORDINATION.md),
which landed as [PR #102](https://github.com/smith6jt-cop/Islet-Explorer-Senior/pull/102).
Keep that handoff's ownership and exchange contract. The October 20 consolidation
plan remains the plan of record; this document records progress and next work.

Evidence checked for this snapshot: fetched Senior refs, the `origin/main` tree,
live GitHub PR metadata, the three open PR heads, and committed handoffs in the
local CODEX checkout. HPG jobs, release payloads and running services were not
inspected in this documentation pass. Run and deployment results below are
attributed to their records, not newly verified here. Recheck before execution.

## Start here

1. Start Senior work from freshly fetched `origin/main` in an isolated branch or
   worktree. At this snapshot it is `1e98a4de48832828be7e2e16c6e7b7e6f8bcf16e`.
   The main workstation checkout is on `maldi/serial-platform` at `00ec032`,
   with untracked `CURRENT_DATA.txt` and `IO60panc2nd Notes.md`. Preserve those
   files and all ignored data; do not reset or clean the checkout.
2. Read root [CLAUDE.md](../CLAUDE.md),
   [typing scope](TYPING_SCOPE_2026-10-26.md), this file and the October 28 handoff.
   A dated narrative, a PR description and a branch head can disagree: inspect
   the relevant code and receipt before treating a capability as available.
3. Senior owns this monorepo and its package, release builder and app. CODEX owns
   its typing run and app. HPG execution owns its jobs and data. Read across
   checkouts as needed, but exchange changes through commits and named artifact
   manifests rather than editing another session's active checkout.

## Merged since the previous handoff

| Change | Evidence and current scope |
| --- | --- |
| Spleen production Protein Measurement | [PR #104](https://github.com/smith6jt-cop/Islet-Explorer-Senior/pull/104), merged as `fbc66a2`. The 12-section primary cohort uses ingest, geometry QC, REDSEA and expression. Its PR records a full HPG run, alignment checks and maintainer gate acceptance on October 31. |
| Spleen per-panel Normalization | [PR #105](https://github.com/smith6jt-cop/Islet-Explorer-Senior/pull/105), merged as `1e98a4d`. PC67's eight sections and PC33's four are fitted separately. Its PR records a completed HPG run and automatic gate acceptance. Unavailable or refused marker fits remain NaN. |

The tracked [spleen primary project](../core/multiplex5000/projects/spleen_primary.yaml)
selects only Protein Measurement and Normalization, with both gates automatic.
The first manual processing stop is Cell Typing. Its intended typing route uses
section-local REDSEA values and levels; normalized units serve region definitions,
classifier rescaling and display. The project does not yet select a production
spleen Cell Typing stage. The old raw-mean pilot is not the primary cohort's
measurement path, and neither measurement nor normalization acceptance establishes
verified cell identity. Spleen remains separate from the first active IO60
pancreas+pLN release scope.

## Open work is not on main

| PR / branch | Observed head and remaining work |
| --- | --- |
| [#103](https://github.com/smith6jt-cop/Islet-Explorer-Senior/pull/103), `redsea/despeckle` | Open, `d292bc2`. Optional median despeckle and retention of no-nucleus cells solely as spillover context; both default off. No-nucleus cells remain excluded from analysis and estimation. Review and merge separately; do not assume either option exists on main or enable it globally. |
| [#106](https://github.com/smith6jt-cop/Islet-Explorer-Senior/pull/106), `maldi/serial-platform` | Open, `00ec032`. Adds the MALDI inventory/ingest and serial registration/integration routes. The PR reports 32 IO60 sections from 18 donors and an inventory gate awaiting tissue/neighbour inspection. Real registration and integration are not established by its synthetic checks. Registration to a Xenium neighbour is unbuilt. |
| [#107](https://github.com/smith6jt-cop/Islet-Explorer-Senior/pull/107), `typing/per-image-method` | Draft, `df3fe1d`. Ports the CODEX per-image method pinned to `578767e` into `multiplex5000.positivity`. Real-data parity remains to be demonstrated; the PR reports that test skipped in its environment. The latest commit adds a PC67 spleen registry with open proposals. |

Two stale descriptions need care when resuming those branches:

- The MALDI plan still says the implementation is uncommitted. It is committed
  and pushed at the head above, but remains unmerged. MALDI relates serial-section
  spots to objects and continuous signed distance, never to the same cells.
- PR #107's body says there is no spleen registry yet. Commit `df3fe1d` adds it,
  including `phenocycler-pc67__spleen__proposals.csv` with P1–P9 still `proposed`.
  `require_reviewed()` refuses open proposals. These are not approved rulings.
  In particular, P7 identifies a needed tissue-reference change: the port's
  hormone-cell checks cannot serve spleen as written. cDC naming, FOXP3 state
  handling and the other proposals also remain review work. No PC33 registry
  is present in this PR.

## CODEX has advanced to v20

The local CODEX checkout at `/data/panc_CODEX/Islet-Explorer` is on
`spatial-review-20261020`, clean at
`307bc1b58afca301e71f1d141dde475b9b6be195`. Its status reports 17 commits ahead
of its existing tracking ref; that ref was not fetched in this pass. The
committed `docs/handoff_20261102.md` and updated `docs/typing_handoff_20261031.md`
report:

- **v20 delivered:** `data/development/uniform-20donor-20261030-v20`, with
  `data/development/qupath-review-20261030-v20`. Final reconciliation passed
  for 20 images and 3,199,606 QC-eligible cells.
- The owner chose DG for Lymphatic objects: despeckled detection plus the PDPN
  level gate. Nerve and Capillary detection retain raw pixels under the
  November 2 ruling. This is not a blanket PhenoCycler despeckle decision.
- The development site serves `pubapps-2026-11-uniform-v20`, app commit
  `c974946`, since November 2 at 06:07 UTC. The record reports reviewer-access
  checks 15/15 and browser checks 49/49 through the UFRC proxy.
- Production has been stopped since October 26; its recorded release is still
  v10. Owner review of v20 and a production decision remain open. Development
  deployment is not production promotion or adoption by Senior.

These records supersede the October 28 snapshot's pending-v18 and served-root
statements. The old `docs/typing_change_notice_20261027.md` still contains those
older claims. The CODEX owner needs to reconcile that notice and publish the
newer committed handoff so remote sessions can retrieve the same evidence.
The November 2 handoff's count of 16 unpushed commits predates its own handoff
commit; local status now reports 17 ahead. Do not modify or push CODEX as part
of a Senior documentation change.

PR #107 still pins method parity to CODEX `578767e`, and its test defaults to
`data/development/uniform-20donor-20261029-v19-qc-work` under the CODEX checkout.
That is a processing run, distinct from the assembled v20 root. Reconcile its
actual provenance with the newer v20 record before running parity. Check the
coverage hashes and selected levels; do not substitute a root or weaken the
provenance assertions merely because v20 exists. A successful parity comparison
would establish transfer fidelity, not biological correctness. The newer
structure/context decisions also need an explicit disposition: PR #107 does
not port construction of context from objects.

## What consolidation still lacks

The [stage registry](../core/multiplex5000/multiplex5000/stages/__init__.py) on
main implements Source Inventory (00), Protein Measurement (03), RNA Ingest
(04), Normalization (05), Typing (06), Registration (07) and Integration (08).
Several are bounded tissue or same-slide routes, not evidence of a completed
cohort-wide release. Image Preparation (01), Segmentation (02), Spatial Analysis
(09), Cohort Analysis (10), Release Build (11) and Release Verification (12)
remain unbuilt there. MALDI Ingest (04b) exists only on PR #106's branch.

The default [IO60 project](../core/multiplex5000/projects/io60.yaml) still
selects only 00, 03 and 05. The dedicated pLN production measurement, full
pancreas+pLN all-modality typing, release build/verification and display-only
app transition remain work to complete and verify. Bounded pilots do not
establish these finish conditions.

For HPG Xenium-only acquisitions 0041323/0041326, the last status available
from the October 28 handoff is accepted Source Inventory, with RNA Ingest and
Typing not yet run. **Current job and artifact status was not rechecked.**
Verify the HPG preservation branch/dirty-source disposition, installed package
commit and scoped project receipts before resuming; do not assume the older
`3809271` wheel pin is still installed.

The research/release distinction in [typing scope](TYPING_SCOPE_2026-10-26.md)
still applies. Separate raw-gene comparisons can proceed, and candidate cell
calls can support clearly provisional comparisons. RNA Ingest acceptance is
required before Typing within the release workflow. Those gates are not a
global hold on research, and two donors cannot establish a population disease
effect. Keep the inclusive gamma-delta screen `(CD3E or CD3G) and (TRDC or TRGC2)`
and report stricter tiers separately. Preserve exact cell keys, nucleus
admission and the recorded unavailable status of damaged 6533 same-slide pLN.

## Next work by owner

1. **Senior typing:** reconcile PR #107's pinned method and processing-run
   provenance with CODEX v20; obtain the manifest/checksum evidence; run parity
   in the matching environment; resolve the spleen proposals and tissue hooks;
   then wire the reviewed method into the intended production stage. Keep the
   PC67 and PC33 panel requirements explicit.
2. **CODEX:** publish the newer handoff and reconcile the stale change notice;
   supply the exact producer, input/output roots, registry, levels, calls,
   exclusions, checks and hashes. Keep its app/data decisions in that session.
3. **Senior/HPG:** verify live Xenium-only receipts and source preservation,
   continue the scoped RNA work, and complete the selected first-release
   pancreas+pLN roster and production routes. Treat the reported spleen runs as
   evidence to inspect, not a reason to repeat completed processing by default.
4. **MALDI owner:** complete the inspection rows, then follow the branch's
   inventory, ingest, registration review and integration sequence. Keep its
   candidate release separate and its Xenium-neighbour limitation explicit.
5. **Consolidation/release:** finish the missing stages needed by the selected
   release, replay into new output folders, verify the candidate and its app,
   and activate through its pinned release pointer. Preserve the two apps and
   separate assay datasets. Retire duplicate code and archive data only after
   provenance and the exact reviewed disposition lists support it.

Validation of this handoff is documentation-only: source/PR cross-checks,
relative-link checks and `git diff --check`. No scientific pipeline, service,
gate or release pointer was changed, and no runtime test result is claimed
from this documentation pass.
