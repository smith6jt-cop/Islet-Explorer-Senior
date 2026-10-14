# HANDOFF 2026-09-18 — same-slide Pheno-Xeno integration (6533, HDL112)

Written at the end of a long session, for a session that opens VSCode in
`/blue/maigan/smith6jt/Islet-Explorer-Senior-main`. Read this before touching anything.

**Nothing in this repo is committed.** 29 files are modified or new. `git status` before you
start. HEAD is `990c5b2`.

---

## 0. READ THIS FIRST — there are TWO forks of the Xenium typing lane

They are both git-tracked, in different repos, and they have diverged:

| | path | size | dated |
|---|---|---|---|
| **monorepo (this repo)** | `analysis/xenium/scripts/hierarchical_phenotype.py` | 64,920 B | 2026-09-15 |
| stale fork | `/blue/maigan/smith6jt/Xenium_Analysis/scripts/hierarchical_phenotype.py` | 62,793 B | 2026-06-12 |

Diverged files and their differing-line counts:
`hierarchical_phenotype.py` 87 · `finalize_panel.py` 40 · `tissue_markers/pancreas.py` 21 ·
`validate_redesign.py` 11 · `state_panels.py` 11 · `tissue_markers/_identity_core.py` 11

The monorepo copy runs in **`scvi-env`**. The `Xenium_Analysis` copy runs in
**`xenium_analysis`**. `analysis/xenium/CLAUDE.md` says the lane was absorbed here 2026-08-28
and that its decisions #1–#17 must not be renumbered.

**I ran the STALE fork for HDL112.** The result below is real and the job exited 0, but it was
produced by the June copy in `xenium_analysis`, not by the September copy in `scvi-env`. Decide
whether to re-run through `analysis/xenium/` before building on it. Reconciling the two forks is
probably the single highest-value cleanup available.

---

## 1. What landed

### HDL112 pancreas — per-cell RNA typing DONE
- `data/raw/0081366.zarr` (5.9 GB) — ingest, SLURM 42503076, 4 min.
- `data/processed/0081366/0081366_phenotyped.h5ad` — SLURM 42554378, **rc=0**, finished
  2026-09-18 09:02. **644,104 cells** (from 671,211 pre-QC) with per-cell `celltype`,
  `celltype_broad`, `celltype_confidence`, `celltype_fine`, `celltype_lineage`.
- Both paths are under `/blue/maigan/smith6jt/Xenium_Analysis/`, NOT this repo.
- Scripts: `Xenium_Analysis/scripts/slurm/112_ingest.sh`, `112_phenotype.sh`,
  `112_phenotype_only.sh`. Logs `Xenium_Analysis/logs/h112_*`.

The first phenotype attempt failed on a **correct fail-loud guard** (decision #12): `AMY1A` is
on the hAtlas custom panel but NOT on the stock Prime 5k panel the Pheno-Xeno sections carry,
leaving `Acinar` with one marker. Fixed by adding `CTRL`, which the vendored identity_core
already marks `subtype_hint=acinar` and which is present on every panel in the project. The
threshold was NOT lowered.

**That fix went into the STALE fork only.** `analysis/xenium/scripts/tissue_markers/pancreas.py:45`
still reads `"Acinar": ["AMY1A", "CUZD1"]` and will hit the same guard on any Prime 5k section.

Also absent from Prime 5k, so the MARKER_DICT thin-fine fallback is thinner there than on
hAtlas: `PTPRC TRAC TRBC2 NKG7 IAPP MAFA NKX6-1 GC CHGA KRT19 ACTA2`. The identity_core itself
is unaffected — 478/478 genes present on both panels.

### Pancreas protein arm — COMPLETE, both donors

`core/phenocycler/data/integration/sameslide/donor_id={112,6533}/roi=panc/serial={0081366,0080870}/`

| | 112 | 6533 |
|---|---|---|
| transform | `refined` PASS 0.0997 µm, 340/359 | `refined_field2` PASS 0.2016 µm, 318/336 |
| joined | 669,610/671,211 with 280 protein features | 712,930/717,461 |
| lifted calls | 659,854, purity 0.982 | 703,641, purity 0.985 |

Chain of custody verified for 112: the sha256 the gate certified is byte-identical to the
`refined.json` on disk, and both `build_joined.json` and `lift_calls.json` record that stem.
`protein_calls.parquet` and `joined_cells.parquet` carry identical `cell_id` sets.

### 6533's transform is now reproducible

`scripts/restain/run_field_refine.py` regenerates `refined_field2` from the repo. Pass 1 matches
the released field to 5.4e-6 µm; pass 2 to 4.6e-5 µm; a fresh gate matches to all 10 printed
decimals (not bit-identical — `median_offset_um` differs in the last ULP,
`median_offset_um_bounded` at the 5th decimal). Three conventions had to be MEASURED: a 125 µm
half-node grid anchor; pass 2 fitted to all 318 `measured` windows (NOT the 317 that also
improve contrast — the contrast filter is what broke every earlier attempt); and pass 2
evaluated at nodes displaced by pass 1 (`q = u + D(u)`), worth up to 1.85 µm.

`refined_field2.*` is sha256-verified untouched. `build_joined.py:57` and `lift_calls.py:61`
still default to it. **Open decision:** the corrected-anchor chain `refined_field4` reaches
median **0.1631 µm** vs 0.2016. Adopting it invalidates every downstream parquet built on
`refined_field2`. Left unadopted deliberately.

---

## 2. THE INTEGRATION PIPELINE HAS NEVER RUN

`phenocycler.integration.pipeline` is the real orchestrator — `ORDER_SAME_SLIDE`, a `sameslide`
stage, `postxen`, `--tissue` defaulting to pancreas **and** pancreatic lymph node together, and
verdict-gating via `qc.ALLOWED_BY_VERDICT`. Two tissues × two modalities on one slide is a
designed feature of it.

```
python -m phenocycler.integration.pipeline --status --donor 6533 --donor 112
```

reports **11 of 15 stages BLOCKED** on one thing, with `cached files ignored`:

```
phenocycler_identity_release = ../../data/development/broadtype/release      MISSING
xenium_identity_release      = ../../data/development/xenium_typing/release  MISSING
```

`config.ini` line above them: *"Optional declared identity releases. Blank is the current
fail-closed state."* `data/development/` contains only `celltype` and `xenium_celltype`.

**Everything produced by hand this session — registration, join, lift, protein calls, decision
records — lives in trees the pipeline explicitly ignores.** `scripts/celltype/run_sameslide_protein.py`
reads `phenotype_workflow.csv` directly and so bypasses the identity-release gate the pipeline
exists to enforce. That is the central unresolved architectural fact.

The maintainer reported (2026-09-17) typing all 21 cohort pancreas donors per-cell —
17,490,305 cells, 15,318,868 called — to `/data/xenium_typing/panc_2026-09-17/release/`. That
path does not resolve on this machine; it is on the box where `/data/xenium` is real. If it is
the `xenium_identity_release`, half the gate is satisfied. **Verify before assuming.**

`section_modes.csv` already declares all four same-slide pairs with evidence: `6533/panc/0080870`,
`112/panc/0081366`, `6533/pln_1/0080870`, `112/pln_1/0081366`.

---

## 3. pLN — state and the one thing that is blocked on the maintainer

**The same-slide pLN Xenium bundles are at
`/orange/maigan/20260910__211837__6533_HDL112_pLNpanc_PhenoXeno/`** (all four: both donors ×
panc/pLN). Not under `Xenium_T1D_PancPLN/`. 6533 pLN = 29,708 cells (1227 × 2007 µm);
112 pLN = 602,113 cells (6942 × 8664 µm).

**6533's pLN is three physical sections**; only one survived to be run. 112's is one intact node.
Measured from the newly uploaded QuPath export: **all 34,841 pLN protein cells are in section 1**
(x 2869–4628, y 1045–2525 µm), zero in sections 2 and 3. So `pln_1` is correct.
`6533_with_pLN.tsv` is a clean superset — 1,062,981 rows = the same 1,028,140 pancreas objects
plus 34,841 pLN.

**BLOCKED ON THE MAINTAINER: 112's pLN protein export does not exist yet.** No
`Cells_pLN_0081366…geojson`; `0112.tsv` is still the 09-15 pancreas-only file (verified: zero
objects below y 8500 µm). Until it lands there is no 112 pLN protein arm.

### 6533 pLN registration — runs end to end, does not pass

Best result: crop taken from where the protein cells actually are
(`--crop-um 2719,895,4778,2675`), median **4.03 µm** against a **1.34 µm** pass bar, 9/14
windows. Output `core/phenocycler/data/development/sameslide_6533/out_pln_s1tight/`.

Ruled out by measurement, so do not re-investigate:
- **Not nucleus detection.** Detection IS poor (`n_if_nuclei` 13,171 vs 34,841 segmented;
  `detection_ceiling` 0.455; `NUCLEUS_UM=6.0`/`MIN_SEPARATION_UM=4.0` are pancreas numbers
  against packed lymphocytes). But re-detecting at 0.5 and 0.25 µm/px changed the result by
  **zero, to 16 decimal places**, because the winning transform is `image_translation`
  (`stages: [valis, whole_section_ncc]`), which is image-based. The nuclei candidates lost to it.
- **Not the crop.** Tightening it bought 0.6 µm and no more.
- **Not instability.** The two held-out halves agree to 0.34 µm on translation.

**Open decision:** a 4.03 µm median with 7.46 µm p90 on a section known to be damaged looks like
the pancreas story (`refined → refined_field → refined_field2`). `run_field_refine.py` exists.
But the pancreas field was fit on **318** measured 500 µm windows and the pLN has **9**, because
the ROI is 2.1 × 1.8 mm. `WINDOW_UM = 500.0` in `run_align.py:51` is a module constant, not a
flag; dropping it to ~150 µm gives ~90 windows but changes what the gate's thresholds mean.
Fit a weakly-constrained field, or shrink the window first? Maintainer's call, not made.

---

## 4. Uncommitted changes (29 files)

### Mine — crop support for a small ROI in a large scan
`scripts/restain/`: `export_dapi.py` (crop + `--crop-um` + per-axis origin),
`registration_inputs.py` (`crop_x0_um…crop_y1_um` columns, all-four-or-none; `as_origin_xy`),
`nuclei.py` (`detect` per-axis; `if_plane` re-applies the crop, else asking for a finer plane
silently swaps the ROI for the whole slide), `refine.py`, `run_image_seed.py`, `run_refine.py`,
`valis_register.py` (resolves the pair inline — it runs in the VALIS venv and imports nothing
from the repo).

**The invariant:** a crop offsets x and y independently, so `export.json`'s `if_origin_um`
becomes a **per-axis pair** when cropped and stays a scalar otherwise. Read it through
`RI.export_origin_xy` / `RI.as_origin_xy`, never `float()`. `if.origin_um` is the effective
origin; `if.origin_um_grid` is the pyramid's own. Uncropped behaviour is unchanged.

`scripts/restain/tests`: **9 failed / 195 passed both before and after** — verified by restoring
the originals and re-running. All nine pre-date this work.

### Mine — bundle resolution
`xenium_celltype/xenio.py` (`load_donor(bundle=)` override) and `run_donor.py` (resolves through
`hpg/sites.py`). `XENIUM_ROOT = /data/xenium` is the other workstation's layout; the nPOD
sections are on `/orange` here, which `hpg/sites.json` already records.

### Agent-authored — reviewed 2026-09-18, sound
One coherent change: **a donor is not a section.** Three independent `.iloc[0]` selections
(`load_donor`, `build_phenotyped_base.build_one`, `hpg/sites.resolve`) all returned the **2024
nPOD run** for 6533 and 112, silently, where the Pheno-Xeno bundle was wanted. Now one fail-loud
`select_bundle_row` with no default. I verified its three self-claims: `serial_of` genuinely
delegates to `xenio`; no stray `.iloc[0]` remains for bundle selection; `h5ad_refusal` is applied
in `load_donor` and surfaced by both consumers (the h5ad carries no serial, so pairing it with a
2026 bundle joins two sections on `cell_id` and would NOT fail, just return nonsense).

New/changed: `hpg/{preflight,sites,stage_call,stage_prep}.py`, `build_phenotyped_base.py`,
`build_decision_records.py`, `decision_record.py`, `pln_background_audit.py`,
`pln_control_sets.py`, `tests/test_bundle_resolution.py`, `tests/test_decision_record.py`,
`run_field_refine.py`, and the CSVs.

**Two defects found in review, both fixed:**
1. The pLN agent added the same-slide pLN rows to `xenium_paths.csv` (correct — both bundles
   verified on disk) but left `test_unambiguous_donors_still_resolve_without_a_serial` asserting
   the old two-pair dual set. Updated to four pairs, still pinned not derived.
2. `hpg/kodama.lock` (untracked, written 2026-09-17 14:35) causes the `test_verify_run` failure.
   With it moved aside: **59 passed, 0 failed**. `verify_run.py:77` reads `HERE/"kodama.lock"` —
   live on-disk state — so check 9 passes or fails depending on whether anyone ran the lane on
   that machine. The lock also records `host` and an absolute `r_libs_user`, i.e. machine state
   inside `scripts/`. **Not fixed** — make the test hermetic, and/or move/gitignore the lock.
   Where it belongs is a maintainer decision.

`scripts/xenium_celltype/tests`: **1 failed / 58 passed** as it stands, 59/0 without the lock.

---

## 5. Known defect NOT fixed — affects the protein calls

`run_step3.build_scale()` computes floor/top for the 7 `panel.NUCLEUS_MEASURED` markers
(FOXP3, Ki67, PCNA, SOX2, TCF_1, TOX, TP63) on the **DAPI ratio**, but
`run_sameslide_protein.py:495` applies that scale to the **raw nucleus means** (medians 78–534
against a 0.09–3.93 scale). After `apply_scale`'s clip, **99.70–99.99% of cells sit at exactly
1.0** — those columns carry no discrimination.

**5 of 37 references depend on them, covering 246,349 of 112's 659,854 lifted calls**, of which
`Exo-acinar-ductal-proliferating` (Ki67 + PCNA) alone is 224,761. This affects the CALLS, not a
figure. A verifier refuted the obvious remedy: dividing by DAPI leaves 5 of 7 markers *below*
their own scale floor.

Related: `Endo-beta-HLA-A-hi` has **0 cells in both donors** and is structurally unreachable
inside islets — `run_islet.subtype_endocrine`'s HORMONE dict is exactly
`{beta: INS, alpha: GCG, delta: SST}`, and the HLA-A-hi variant is in the 11-reference set
`intruder_references` excludes from step 4 of the core cascade. The T1D MHC-I phenotype is not
being measured.

---

## 6. Corrections — things asserted wrongly this session

Recorded so they are not inherited.

1. **"Your repo types Xenium cells per-cell"** — over-generalised from four NON-cohort bundles
   (0041323, 0041326, 0041435, 0041534) in `Xenium_Analysis/data/processed/`. No cohort donor
   had a per-cell RNA lineage before 2026-09-17. `run_donor.py` is honest about its scope:
   `DEFAULT_MAX_CELLS = 50_000`, *"The deliverable is donor-level COMPOSITION, not a label on
   every cell."* The gap was real, not self-inflicted by lane choice.
2. **"No same-slide pLN Xenium bundle exists"** — wrong. Searched `Xenium_T1D_PancPLN/` and
   `Xenium_Analysis/data/` and stopped. They are at `/orange/maigan/20260910__211837__…`.
3. **"1 failed / 58 passed is a pre-existing baseline"** — true relative to the agents' edits,
   false relative to a clean checkout. See §4 defect 2.
4. **6533's `refined` stem FAILs at 5.47 µm**, not 38.798 µm; that figure is a superseded
   2026-09-15 file in the second gate tree.
5. `joined_cells.parquet`'s 280 protein features are **4 compartments × 70 markers**
   (Cell/Cytoplasm/Membrane/Nucleus), not 280 markers.
6. Terminology: the 70-plex PhenoCycler DAPI channel is **not** "IF". In this repo `IF` /
   `SLIDE_XENIUM` means the post-Xenium INS/GCG restain — a different image of a different
   section (`phenocycler/sections.py`). Say "PhenoCycler DAPI".

---

## 7. Environment

- Python: `/blue/maigan/smith6jt/miniforge3/envs/xenium_analysis/bin/python3` by **absolute
  path**. `analysis/xenium/` wants **`scvi-env`** instead.
- R: `export PATH="$PATH:/apps/R/4.6/bin"` — **APPEND, never prepend**; prepending shadows
  conda's `python3` with an R-bundled one that has no pandas. `R_LIBS_USER=/blue/maigan/smith6jt/R_libs`
  (KODAMA 3.3, Rtsne, umap installed 2026-09-17).
- VALIS: `~/.venvs/valis/bin/python` (wrapper sets its `LD_PRELOAD`).
- SLURM: GPU work goes to `--partition=hpg-b200 --gres=gpu:1 --qos=maigan`. This session ran on
  a CPU node with no GPU visible.
- Data: `/blue/maigan/smith6jt/Xenium_Analysis/data`. A path written `/data/...` means THAT,
  never the filesystem root.

**Hazards that cost time this session:**
- NEVER `pgrep -f` / `pkill -f` a pattern contained in your own command line — it matches itself.
  Use `ps -eo pid,cmd | awk '!/awk/ && /pat/'` and kill explicit PIDs.
- `cmd | tail` returns tail's exit status. Capture `${PIPESTATUS[0]}`.
- `cd X && cmd &` backgrounds the whole list, so the parent's cwd never changes.
- The QuPath TSVs are 5–7 GB. Use `awk`, not pandas.

---

## 8. Suggested first moves

1. `git status` — 29 uncommitted files. Nothing is pushed.
2. Decide the **fork question** (§0). Everything else is downstream of it.
3. Verify whether `/data/xenium_typing/panc_2026-09-17/release/` satisfies
   `xenium_identity_release`. If it does, the integration pipeline may be one release away from
   running for the first time.
4. Port the `Acinar` + `CTRL` fix into `analysis/xenium/scripts/tissue_markers/pancreas.py:45`
   if the monorepo lane is canonical.
5. The two open decisions: adopt `refined_field4` (§1), and the pLN field-vs-window question (§3).
6. The maintainer proposed a genuinely useful check: run `hierarchical_phenotype.py` on cohort
   donors as an **independent second method** against the 2026-09-17 release, the role the
   scANVI comparison played (0.814/0.832 on 6436/6442). Feasibility verified — all 478
   identity_core pancreas genes are present on the Prime 5k panel, and the cohort h5ads carry
   `cell_area`, `donor_sex`, `segmentation_method`. It needs `layers['counts']` +
   `layers['lognorm']` (`hierarchical_phenotype.py:1093` exits without them) and
   `var['panel_for_embedding']` from `finalize_panel.py` (line 261).
