# Handoff: Senior monorepo and CODEX Islet-Explorer

Snapshot: 2026-10-28. This is a coordination record for the Claude agent working
in this repository and the Claude agent working in
`/data/panc_CODEX/Islet-Explorer`. Recheck live status before changing either
checkout. The 2026-10-20 consolidation plan remains the plan of record; this
file gives the next work and the exchange between agents.

## Start here

1. **Senior agent:** work from current `origin/main` in an isolated branch or
   worktree. The workstation's `main` checkout was 51 commits behind at this
   snapshot and holds untracked `CURRENT_DATA.txt` and `IO60panc2nd Notes.md`.
   Do not reset, clean, or overwrite that checkout. Read root `CLAUDE.md`,
   `docs/TYPING_SCOPE_2026-10-26.md`, and this file.
2. **CODEX agent:** continue its approved final typing run using
   `docs/final_typing_handoff_20261028.md`, `docs/typing_rulings_ledger.md`,
   `docs/final_typing_plan_20261028.txt`, and
   `output/typing_final_checks_20261028/README.md` in its own repository. Its
   current branch is `spatial-review-20261020`; it has many uncommitted research
   files. Preserve its work before any checkout change. Do not make the Senior
   agent infer readiness from a commit title or an older output root.
3. **HPG execution, if a separate agent owns it:** move new code work to
   `/blue/maigan/smith6jt/Islet-Explorer-Senior-main`. Keep the raw Xenium and
   PhenoCycler payloads at their recorded `/orange` and `/blue` paths, supplied
   through the site file. First preserve the one unpushed commit and dirty code
   in `/blue/maigan/smith6jt/Xenium_Analysis` on a pushed preservation branch.
   Its existing processed outputs are audit context, not inputs to the new
   typing producer. Do not move data or overwrite an active checkout to make
   the code switch.

## Current facts to reconcile, not assumptions of completion

| Area | Verified at this snapshot | Next owner |
| --- | --- | --- |
| Senior code | `main` at `40fde07` (PRs #100 and #101 merged; 19 checks each passed). No open PRs. Five remote branches remain: `main`, two active work branches, two preservation branches. | Senior |
| Stage implementation | Built: 00, 03, 04, 05, 06, 07, 08. Stage 01/02 and 09–12 are not built. The default `io60` project still names only 00, 03, 05; bounded typing projects do not constitute a full active release. | Senior |
| HPG Xenium-only 0041323/0041326 | Project `io60_xenium_only_6512_6618`: Stage 00 complete and accepted; Stage 04 and 06 not run. The HPG installed wheel is pinned to `3809271`; the Senior HPG checkout had reached `main` but has an untracked γδ/islet handoff. | Senior/HPG |
| CODEX typing | Uniform 20-donor v18 has **not** been declared assembled or ready. The live agent has committed canary/identity/acceptance-check code (`c5459f4`) and is refining it. `docs/typing_change_notice_20261027.md` still says not to adopt a root, levels or labels until v18 is assembled and the notice changes. The served app remains on an older root. | CODEX |
| Senior branches | PR #39 was a stale 2026-10-21 execution snapshot; it is closed and its exact tip is `archive/pr-39-20261028`. Merged branch tips were tagged before deletion. Active and preservation branches were retained. | Senior |

The HPG Stage 00 acceptance is signed `Codex` and covers source inventory
only. It is neither typing acceptance nor biological approval. Its accepted
receipt is available for the scoped release workflow; do not turn a question
about its signer into a hold on separate research analysis. The HPG γδ handoff's
suggestion to wait for Stage 06 acceptance before running the final comparison
predates PRs #100/#101 and is superseded for provisional research work.

## One authority and two assay profiles

MultiPlex5000 in `core/multiplex5000` is the maintained processing code. It
will serve the PhenoCycler and CODEX differences through explicit profiles,
not through a copied second typing implementation. CODEX's app and Senior's
app remain separate display clients. The legacy CODEX and PhenoCycler
measurements remain separate datasets and must not be pooled into one cell
table or treated as one assay. One maintained package does not mean one
cutoff for every tissue or assay.

The Senior agent owns changes to this monorepo, its package, project files,
release builder and Senior app. The CODEX agent owns its current final typing
run and CODEX app until it publishes the handoff below. Each agent may read the
other's work, but neither edits the other's active checkout or data root.
Exchange through a committed notice and a named immutable output root, not a
moving `current` path. If HPG execution is a third session, it owns only HPG
jobs/data and submits code changes to the Senior monorepo.

## Immediate work, in parallel

### CODEX agent: finish and publish the typing handoff

1. Complete the already approved CODEX final typing plan. The owner chose the
   per-image SMA negative-population p99 and CD20 epithelial-negative p99.9
   options recorded in `docs/final_typing_handoff_20261028.md`. Finish its
   Step 4 levels, Step 5 canary, Step 6 identity/tiles and Step 7 report in
   the prescribed order. The existing handoff says to stop for owner review
   at Step 7 before assembling a root or touching QuPath/app deployment.
2. Reconcile the newer code and results with the older
   `docs/typing_change_notice_20261027.md`; update the notice only when the
   final root is actually assembled. Publish a manifest naming the exact
   CODEX input root, final output root, producer commit, panel/marker registry,
   positivity levels, per-cell calls, label dictionary, controls, exclusions,
   counts, checks and SHA-256 values. Record which rules are CODEX-only.
3. Give the Senior agent the committed manifest and evidence report. Keep the
   CODEX app pinned to its current release until the final result is reviewed
   and a separate display switch is ready. Do not copy a provisional label
   root into Senior or edit Senior's PhenoCycler cells.

### Senior agent: continue without waiting for CODEX v18

1. Preserve HPG `Xenium_Analysis`'s unpushed/dirty source, then move new HPG
   code work to Senior. The HPG checkout can reference existing external data
   paths. Its legacy notebooks and 0041323/0041326 processed outputs remain
   audit evidence; port only unique, reviewed code into Senior. The HPG
   γδ/islet plan is currently an **untracked** file at
   `docs/HANDOFF_2026-10-28_XENIUM_ONLY_GD_ISLET.md` in the HPG Senior checkout;
   preserve/review it before relying on it as a shared instruction.
2. For the two Xenium-only pancreas sections, run Stage 04 and Stage 06 in
   their scoped project when useful for class-specific comparisons. Stage 04
   gate acceptance is a dependency *within that release workflow*. Separate
   research comparisons of available raw counts can proceed now; Stage 06
   candidate calls can be analyzed before its gate is accepted if labeled
   provisional. Do not wait for protein typing, pLN pairing or stages 09/10.
   Those analysis stages are not built. PRs #100/#101 and
   `docs/TYPING_SCOPE_2026-10-26.md` clarify the distinction.
3. For γδ screening, retain the owner's inclusive rule: (CD3E or CD3G) and
   (TRDC or TRGC2), including cells expressing both alternatives. Report
   stricter TCR-chain and background-adjusted tiers separately. Do not limit
   the screen to the old `Immune_T` label or claim a productive receptor from
   Xenium counts. Use exact cell keys, exclude no-nucleus cells, and estimate
   background per section. With one Aab+ donor and one T1D donor, report
   observed contrasts and within-section nulls, not a population disease
   effect. Do not repeat the legacy islet-distance bug.
4. Independently advance package/release work: declare the full first-release
   pancreas+pLN acquisition scope, finish missing Stage 01/02 and 09–12
   functions as needed, replay selected producer outputs in new folders,
   complete the pilot and full rebuild, verify receipts and manifests, and
   prepare the display-only Senior app against a pinned candidate. The bounded
   Stage 06 pilots are evidence, not an active release. Keep 6533's damaged
   same-slide pLN analysis marked unavailable and nonblocking. Spleen remains
   a first *pilot* extension; it does not expand the first active release or
   authorize changes to FollicleCounterQP.

## The CODEX-to-Senior exchange

Once the CODEX agent's final notice says the result is ready, the Senior agent
opens a focused PR that ports the shared **method and rulings** into
MultiPlex5000 with CODEX-specific profile settings. Already recorded Senior
rules include whole-cell fallback when a cytoplasm measurement is undefined,
the CD8/CD4 same-unit comparison and full-profile tie review, and CD68/CD163
macrophage profile labels. The new SMA, CD20, canary and marker-level decisions
need assay-by-assay review before generalizing to PhenoCycler or pLN. Missing
measurements remain unavailable, and a cell without a usable nucleus remains
ineligible. `cell_class` and `refined_class` are successive steps of one
hierarchy, not competing releases.

The PR must name the source CODEX commit and manifest, map marker roles to
`core/panels`, distinguish profile choices from shared logic, and replay a
small fixed sample and the selected output against the published receipt.
Matching old code proves transfer fidelity, not biological correctness. The
CODEX agent reviews the mapping; the Senior agent owns the PR and tests.
After the shared code is in Senior and the CODEX app reads a pinned release,
retire the copied processing code in CODEX on a separately reviewed list.

## Finish conditions for consolidation

- A neutral-directory install of `multiplex5000` and `mplx5k` runs the same
  scientific stages locally and on HPG Slurm; data roots come from site
  settings, not install location or checkout paths.
- One selected pancreas+pLN candidate passes source, panel, identity, count,
  geometry, typing-control and app-render checks and is made active through
  its single release pointer. Every served result names its release and
  evidence limits. CODEX data have their own pinned profile/release and are
  never merged into PhenoCycler measurements.
- Both apps are display-only and read pinned releases. Archive/delete old
  code and derived payloads only after each unique source is accounted for,
  current receipts no longer refer to archive/scratch, and the exact deletion
  lists are reviewed. Raw acquisitions remain referenced by identity/digest.

Do not create a new user-facing hold for ordinary analysis, tests, or porting.
Document remaining uncertainty in the output. Release-backed claims and
activation still use the release checks; a passing CI run or a preliminary
cell label alone is not that evidence.
