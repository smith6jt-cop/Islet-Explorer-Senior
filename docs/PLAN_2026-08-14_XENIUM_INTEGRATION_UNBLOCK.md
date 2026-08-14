# Plan — unblocking PhenoCycler ↔ Xenium integration

**Status:** approved in principle, **not started**. Held pending completion of the
Phase 2 monorepo consolidation (session `monorepo-retirement-commit`).
**Date:** 2026-08-14
**Decision owner:** maintainer (smith6jt)

---

## Decision

Do **not** wait for `core/phenocycler` pair-based typing. Supply PhenoCycler cell
identities to the integration layer from the **existing phenotyped AnnData**, joined by
QuPath UUID, as an explicitly-stamped legacy source.

This is a deliberate departure from `export_pheno.py`'s fail-closed design. It is
compatible with the `CLAUDE.md` non-negotiables **only** under the constraints in
[Hard constraints](#hard-constraints) — in particular, legacy labels flow *downstream into
integration only*, and never back into the source workflow.

---

## The constraint that decides the scope of this work

The legacy phenotypes cover **five** of the twenty run donors, and **none of them are ND**.

| set | n | donors |
|---|---|---|
| current run (`3bf4fc8e391d007766d1`) | 20 | 115, 6374, 6380, 6436, 6442, 6450, 6461, 6476, 6482, 6493, 6505, 6521, 6523, 6534, 6539, 6551, 6566, 6591, 6593, 6623 |
| `xenium_paths.csv` | 26 | above + 112, 6347, 6457, 6533, 6579, 6580 |
| phenotyped H5AD | 15 | 112, 6356, 6450, 6479, 6505, 6516, 6521, 6533, 6534, 6538, 6548, 6549, 6550, 6551, 6563 |
| **all three** | **5** | **6450 (Aab+), 6505 (Aab+), 6521 (Aab+), 6534 (T1D), 6551 (T1D)** |

Every phenotyped H5AD on disk is the same 15-donor Junior cohort
(`panc_CODEX.h5ad`, `CODEX_panc_scvi_2025Nov12.h5ad`,
`CODEX_scvi_BioCov_phenotyped_newDuctal*.h5ad`, `islet_explorer.h5ad`). There is no
larger one. So five donors is a **ceiling, not a starting point**.

**Consequence for the study question.** The question is the immune environment around
islets across ND → Aab+ → T1D. The identity-dependent lane would have 3 Aab+ and 2 T1D
donors and **no ND baseline at all** — no progression axis, and no way to distinguish a
disease effect from the batch confound already recorded in `gate5-batch-confound`.

**Therefore:** identity-dependent cross-modality claims are *not* the deliverable of this
plan. They are a five-donor methods demonstration. The deliverable is the
**identity-independent registration lane across all 20 donors**, which the legacy labels
are not needed for and which is a prerequisite for everything that comes after verified
typing. The legacy-identity lane rides along on 5 donors as a wiring test.

If the ND gap is unacceptable, the alternative is to phenotype the remaining 15 run donors
by whatever method — which is the pair-based typing work this plan was meant to avoid.
That tradeoff is the maintainer's call and is flagged, not decided, here.

---

## Resolved by investigation (2026-08-14)

1. **Which 20 of the 26 donors are the run cohort — answered.** All 20 run donors appear in
   `xenium_paths.csv`; 41 of its 70 rows cover them. The cohort is no longer unrecorded.
2. **The join key exists.** The phenotyped H5AD carries `Object ID` holding the QuPath UUID
   (e.g. `a6f4eacd-f4af-4d98-8154-648b733f6131`), matching `object_id` in
   `runs/<id>/expression/`. **No centroid matching or donor-name parsing is required or
   permitted.**
3. **`expression` already carries everything the geometry lane needs:**
   `donor_id`, `object_id`, `image`, `X_centroid`, `Y_centroid`, `cell_region`, `islet_num`,
   `qc_analysis_eligible`, `qc_estimation_eligible`, and one authoritative intensity for each
   of 68 markers. 1,030,803 cells for donor 115 alone.
4. **`contract.py` over-requires identity.** `REQUIRED_COLUMNS` includes `lineage_native`,
   and `pipeline.py` blocks 13 stages under one flat `PHENO_DEPENDENT_STAGES` set. Stages
   needing only coordinates are blocked alongside stages needing calls. That conflation is
   the artificial half of the block.

## Open, unresolved

5. **Xenium root prefix — mostly addressed 2026-08-14; confirm on HiPerGator.**

   *Correction to this plan's first draft: the claim that "14 rows have an empty
   `source_path`" was wrong — an artifact of how I bucketed the paths, not a property of
   the file. Those 14 were all-blank spreadsheet spacer rows carrying no donor and no roi.
   The manifest was ~95% complete for the cohort all along, and the case for deleting it
   rested on that bad number.*

   Actual state: 56 real rows, 41 covering the 20 run donors, of which exactly **one** has
   no path — donor **6591 `pln_1`**.

   Done: spacer rows stripped; the 6591 gap now surfaces as a warning and a `no_path`
   status instead of being silently filtered; `xenium_root` set; `--verify-paths` added so
   the file is checked rather than merely trusted.

   Still open: `xenium_root = /orange/maigan/Xenium_T1D_PancPLN` is **provisional**.
   `_rewrite_root` substitutes everything up to and including the literal `/xenium/`
   segment, so this value assumes the `HuPrime5kPanTissue/<batch>/<run>/` levels survive
   beneath the new root. Not verifiable here — `/orange` is not mounted. Settle it with
   `python -m phenocycler.integration.manifest --verify-paths` on HiPerGator, then fill or
   formally retire 6591 `pln_1`.
6. **Where this executes.** `/orange` is HiPerGator-only. Either run there, or transfer.
   `/data/QP_Xenium` holds only TMA172/HDL bundles — unrelated samples, not this cohort.
7. **`donor_overrides.csv` is still empty** — `0041323__panc` and `0041326__panc` have
   finished h5ads and no donor.
8. **Legacy phenotype vocabulary needs a crosswalk.** 21 categories
   (Acinar, Macrophage, Structural, Blood Vessel, Ductal, ECAD+, Endothelial, APCs,
   Lymphatic, Alpha cell, Neural, Beta cell, SMA+, Endocrine, CD8a Tcell, **Unknown**,
   T cell, B cell, Immune, Delta cell, CD4 Tcell) against the 8-class common vocabulary in
   `vocab.py`. `Unknown` (8,821 cells) must map to `Unavailable`, never to a lineage and
   never dropped silently.

---

## Architecture

Split the single fail-closed gate into two, along the line that actually matters:
**what a stage reads**, not which workflow produced it.

```
                    ┌─ geometry + intensity ─────────────┐
runs/<id>/expression│  object_id, x/y, area, cell_region,│──► contract (no identity)
                    │  islet_num, 68 marker intensities  │      lineage_native = NULL
                    └────────────────────────────────────┘
                                                              ▲
legacy H5AD ──join on Object ID == object_id──► identity ─────┘
                                                lineage_source = 'legacy_h5ad'
                                                lineage_provenance = <path>@<sha256>
```

`export_pheno` becomes two functions behind one stage:
`export_geometry(run_id)` — always available — and `attach_legacy_identity(...)` — available
only when explicitly configured, and always stamping provenance.

### Stage disposition

| stage | reads | disposition |
|---|---|---|
| `manifest` | CSVs | already runs |
| `import_xenium` | Xenium bundles | already runs |
| `export_pheno` | `runs/<id>/expression` | **rewrite** — geometry always; identity opt-in |
| `structures` | `islet_num` / `cell_region` | **unblock** — QuPath geometry (see below) |
| `register` | coordinates, rasters | **unblock** — all 20 donors |
| `transform` | transforms | **unblock** — all 20 donors |
| `match` | islet polygons both sides | **unblock** — all 20 donors |
| `qc` (registration gates) | Dice, islet RMSE | **unblock** — all 20 donors |
| `figures` (overlays) | registration | **unblock**; composition figures stay identity-gated |
| `vocab` | lineages | identity-gated — marker↔gene half can split out (panel metadata, no calls) |
| `donor` | lineage fractions | identity-gated — **5 donors, no ND** |
| `grid` (niches) | kNN composition class | identity-gated — 5 donors |
| `cellmatch` | endocrine identity | identity-gated — 5 donors |
| `crossmodal` | anchor identity | identity-gated — 5 donors; INFERENCE, QC-gated |
| `postxen` | INS/GCG positivity | **stays blocked** — needs a *call*, not an intensity. The re-stain is a real measurement, but scoring it requires a hormone threshold nobody has verified |

### The `structures` question

Point-set registration (Track A) needs PhenoCycler islets. The current stage seeds DBSCAN on
hormone calls, which do not exist. QuPath already supplies `islet_num` and `cell_region`
per cell.

**Use them, as registration geometry only.** `CLAUDE.md` bars QuPath annotations from being
*training, tuning, or release-gating* inputs; registration is none of the three. Record the
source in the transform's `qc.json`. Never propagate `islet_num` as identity.

The fallback, if this is judged too close to the line, is image-track-only registration
(DAPI↔DAPI / density rasters). That costs Track A, its cross-check against Track B, and
makes the islet-RMSE gate undefined — leaving tissue Dice as the sole gate. Not recommended.

---

## Hard constraints

1. **Never write legacy labels into `core/phenocycler/data/runs/`.** Runs are immutable and
   content-addressed. Integration output is a separate root.
2. **Every legacy-derived column carries provenance** — `lineage_source='legacy_h5ad'` plus
   the source path and content hash — in the parquet and in every figure caption and report
   built from it. A reader must never have to ask whether a label is current.
3. **Legacy labels never flow back upstream.** They do not define thresholds, pairs, or
   gates, and do not tune any source-workflow parameter. The prohibition in `CLAUDE.md`
   stands unchanged.
4. **Join on `object_id` only.** Never centroid matching, never donor-name parsing.
   Fail loudly on an unmatched UUID rather than falling back.
5. **`Unknown` → `Unavailable` with a recorded reason.** No silent drop, no "Unassigned"
   bucket.
6. **A missing measurement is unavailable, not zero.**
7. **Donors 6457 and 6579 stay excluded** (`cohort.py`). Both appear in `xenium_paths.csv`;
   the manifest must drop them.
8. **Any identity-dependent result is labelled `n=5, no ND` wherever it appears.**
9. **No result is reported without the figure that shows it** — per
   `provide-data-not-just-assessments`.

---

## Work plan

**Phase 0 — wait.** Phase 2 consolidation finishes and lands. Integration code is at
`core/phenocycler/phenocycler/integration/`. Do not edit before then.

**Phase 1 — data reachability.** Resolve item 5: reconcile `/orange/maigan/Xenium_T1D_PancPLN`
against `xenium_paths.csv`, fill the 14 empty `source_path` rows, decide HiPerGator vs local,
fill `donor_overrides.csv`. Then `manifest --report` should show 20 paired donors and zero
`donor_unknown` among them. *No code changes — this is a data-truth phase, and everything
downstream is worthless if the manifest is wrong.*

**Phase 2 — the geometry split.** Make `lineage_native` optional in `contract.py`; split
`PHENO_DEPENDENT_STAGES` into `PHENO_IDENTITY_STAGES` and geometry stages; write
`export_geometry` against `runs/<id>/expression`; extend the test suite (currently 200 tests,
all synthetic) to cover a geometry-only contract table. Gate identity stages on the actual
presence of `lineage_native`, not on a hardcoded list.

**Phase 3 — registration, all 20 donors.** `structures` → `register` → `transform` → `qc`.
This is the expensive, risky, identity-independent core, and the real prize: when verified
typing eventually lands, only the identity stages remain. Deliverable is the registration QC
table and per-donor overlay figures, reviewed against the measured floor per
`review-figures-must-show-background-honestly`.

**Phase 4 — the legacy identity lane, 5 donors.** `attach_legacy_identity`, the 21→8
vocabulary crosswalk, then `vocab` / `donor` / `grid` / `cellmatch`. Treated and labelled as a
wiring test, not a cohort result.

**Phase 5 — decide.** With registration proven on 20 donors and the identity lane proven on 5,
the question becomes concrete: phenotype the other 15 run donors by some verified method, or
accept a 5-donor identity result. Revisit then.

---

## Notes

- `postxen` (S1c) is the one place protein and RNA are measured on the same cells. It stays
  blocked here, but it is the highest-value thing a verified hormone call would unlock —
  worth remembering when pair design resumes.
- Two existing runs (`3bf4fc8e391d007766d1`, `c8412757d12de36ca056`) contain populated
  `assignments/` and `marker_evidence/` from the retired `ordered_residual_typing` lane.
  These are **not** an identity shortcut. If a future reader is tempted, `run.json` stamps
  them `production-default-restore-ordered-residual-typing-v1`.
