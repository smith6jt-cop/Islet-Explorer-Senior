# Handoff — Xenium data, cell typing, and the integration rework

> **SUPERSEDED IN PART — read this first, 2026-09-18.** This handoff argues an acinar/ductal axis and states "The Exocrine core is already an acinar call driven by CUZD1, not a ductal one" (§2b, and again around the CUZD1/TFF1 discussion). **That claim is not supported and is not the established position.** The Xenium 5K panel carries no PRSS1, CPA1, CELA3A, CTRB1, CTRB2, AMY2A or REG1A — measured per section — so **Exocrine is ONE label and no acinar and no ductal claim is made or supported**; see `docs/PANC_XENIUM_TYPING_2026-09-17.md` §5 and the same rule for the 60-plex in the root `CLAUDE.md`. Sections written as live 2026-08-27 CORRECTIONS in this file read as current and are not. **Two
> further places carry the same retracted axis as live work and are covered by this banner:** §4's
> "Do not report these numbers" reason 1 ("The acinar no-axis problem (§2b) is unaddressed ... most
> likely by `Exocrine` (the nearest, ductal)") and §6's next step 2 ("Solve the acinar no-axis
> problem"). There is no acinar axis to solve for and no ductal reading of Exocrine; the pilot
> numbers stay unreportable for the reasons in §2a and §2e, which are untouched.

**Date:** 2026-08-26
**Branch:** `feat/dct-step0-pinning` (HEAD `c9edcf0` at time of writing; the peer
`marker-qc-uniform-preprocessing` session has been committing to the same branch throughout)
**Plan file:** `~/.claude/plans/wobbly-wibbling-quilt.md`
**Nothing in this session is committed.** See [Uncommitted state](#uncommitted-state).
*(2026-09-02: since committed as `d0b756f` … `3f68e24`; the branch is 56 commits ahead of `main`.
Current state: `docs/WORKFLOW_2026-09-02_MAP.md` § 2.4.)*

---

## 1. What landed

### 1a. `/orange/brusko` purged from the repo — DONE

`/orange/brusko/xenium/...` is another lab's copy. Removed everywhere; **0 references remain**.

- `data/integration/xenium_paths.csv` — 55 rows rewritten to **site-neutral relative paths**
  (`HuPrime5kPanTissue/<batch>/<run>/output-XETG...`). Verified by `git diff` that every
  changed line differs *only* by the stripped prefix, so donor ↔ serial ↔ run is untouched.
  This is safe because `parse_bundle_path` (`manifest.py:146`) matches with `.search()`
  anchored at the end of the string — the prefix contributes nothing to identity.
- `manifest._rewrite_root` — now resolves relative paths against `xenium_root`, keeps the
  `/xenium/` splice for un-migrated absolute rows, and adds `_find_bundle`, a bundle-name
  search that **refuses on an ambiguous match** rather than guessing.
- `verify_bundle_paths` — two new statuses: `no_root` (relative paths, no root configured — a
  config gap, not a missing file) and `relocated` (found by name, not where recorded — data
  fine, manifest stale).
- `config.ini` — `xenium_root = /data/xenium`; comment rewritten to record that data arrives
  by Globus **transfer**, not mount.
- 4 new tests. Full suite **574 passed, 2 skipped**.

### 1b. The record corrected — DONE

Several documents generalised a single model's failure into a claim about the platform. Each
is now scoped to what was measured. **Nothing measured was removed; only the reach of the
conclusion changed.**

- `docs/DEEPCELL_TYPES_REVIEW_2026-08-24.md` — "The immune compartment is not delivered" →
  "This model, at this checkpoint, did not recover the immune anchors", plus a dated standing
  note explaining the rescope and closing with *do not cite this document as a capability
  verdict on the platform or the panel*.
- `scripts/celltype/README.md` — the CD11c/CD209 line restated as a measurement about two
  markers' bright tails, with the mechanism, not a verdict on the panel.
- `METHOD_HISTORY.md` needed no change (already scoped to gate order).
- Memory: `gate5-batch-confound` extended with the crossed-batch addendum; new
  `no-capability-verdicts-from-one-result` rule written and indexed.

### 1c. Xenium data transferred — DONE

Globus Connect Personal (already registered here) → `/data/xenium/`. **Both transfers
SUCCEEDED, 0 faults**, independently verified against the remote listing (no missing files, no
extras, no size mismatches).

| stage | files | size | what |
|---|---|---|---|
| 1 | 46 | 50.8 GB | `h5ad_out/` — all per-donor AnnData |
| 2 | 60 | 7.7 GB | `cells.parquet` + `cell_feature_matrix.h5` + `cell_boundaries.parquet`, all 20 run-cohort pancreas bundles |

`manifest --verify-paths` now reports **`ok` for exactly the 20-donor run cohort**, all
`roi=panc`. A bundle is ~14 GB but only ~385 MB of each was taken; `morphology.ome.tif`
(10.3 GB) and `transcripts` (3.3 GB) are not needed unless the image registration track
underperforms.

**The four donors missing from maigan were in brusko**, as the maintainer said: MM046
(6450, 6505) and MM047 (6534, 6623) — 2 of the cohort's 3 Aab+ donors. They were pulled from
there; the other 16 from maigan. **No brusko path entered the repo or config** — the local
copy under `/data/xenium` is ours, paths are relative, and reading brusko once to make that
copy is acquisition, not a recorded dependency.

---

## 2. What was measured, and what it refutes

### 2a. The per-cell identity-core argmax is refuted — do not ship it

Donor 6539, 50,000 cells, against random gene sets matched on core size **and per-gene
detection decile**, 10 draws:

| lineage | genes | REAL % | null mean | null p5–p95 | |
|---|---|---|---|---|---|
| Endocrine | 18 | 2.83 | 5.84 | 4.49–7.66 | **below null** |
| Exocrine | 8 | 75.46 | 61.18 | 38.68–77.76 | inside |
| Immune_B | 58 | 1.91 | 6.74 | 2.65–12.87 | **below null** |
| Immune_T | 129 | 1.08 | 5.09 | 0.75–18.50 | inside |
| Myeloid | 111 | 2.78 | 3.08 | 1.53–4.56 | inside |
| Neural | 17 | 1.23 | 1.42 | 0.95–1.92 | inside |
| Stromal | 48 | 7.02 | 4.36 | 2.47–6.37 | above |
| Vascular | 87 | 7.68 | 12.30 | 5.51–34.87 | inside |

**Five of eight are indistinguishable from random gene sets; two do worse than random.** At
~96% per-cell zeros the argmax reports which core holds a high expresser — a property matched
random sets share. Independently, an adversarial review found a size-fair z-score gate still
calls 6.8–13% of cells on pure noise cores.

### 2b. ~~Acinar has no axis~~ — CORRECTED 2026-08-27

> **This section was wrong and is retained only so the correction is legible.** It claimed
> `02_pancreas_exocrine` is "ductal" and that acinar therefore has no axis, and called that
> "the single most important open problem". The panel file contradicts it: of the eight genes
> `CFTR, CLDN1, CTRL, CUZD1, HNF1B, MUC1, SOX9, VTCN1`, **`CTRL` and `CUZD1` are declared
> `subtype_hint=acinar`**, and both are on the panel.
>
> Measured on donor 6539, all 317,610 cells: CUZD1 detected in **65.8%** (mean count 10.69,
> mean lognorm 1.494), CTRL in 40.7%; `CUZD1|CTRL` covers **66.5%** of cells, and the two
> supply **71%** of the core's mean signal. The Exocrine core is already an acinar call driven
> by CUZD1, not a ductal one.
>
> The narrow, true claim: the canonical **zymogen** set is absent (PRSS1/PRSS2/CPA1/CPB1/
> CELA2A/CELA3A/CTRB1/CTRB2/PNLIP/CLPS/AMY2A/REG1A/GP2/PTF1A/RBPJL), as are INS/GCG/SST/PPY.
> Acinar identity here rests on two genes whose acinar specificity **in this data** has not
> been measured — a measurement to make, not a capability to deny
> (see the `no-capability-verdicts-from-one-result` rule).
>
> The real distortion of the pilot shares is §2e below, not this.

### 2c. Zero-evidence cells (fixed, but know the trap)

0.21%–1.23% of cells per donor carry **no identity-core transcript at all**. A naive `idxmax`
assigns every one to `Endocrine` — because it sorts first — fabricating calls in exactly the
compartment islet work depends on. Plus 0.05%–0.14% exact ties.

### 2d. The two batch structures are crossed — the reason for all of this

Every ND donor is in PhenoCycler batch 1 (batch 2: 2 AAB, 5 T1D, **0 ND**). But **6 of the 11
Xenium instrument runs contain both an ND and a non-ND donor**, and 3 span both PhenoCycler
batches. An effect present in both modalities cannot be explained by either instrument's batch
structure. This needs no registration and no cell matching, and survives QC verdict FAIL.

---

### 2e. What actually distorts the pilot shares (measured 2026-08-27)

Three defects, all verified against `group_calls.csv` and donor 6539's own matrix. Together
they explain the k-instability that §2b was blamed for.

**`Immune_B` is a UPR readout, not a B-cell measurement.** `04_immune_B_plasma` is the only
identity core with substantial state-module contamination. Its 9 ER/secretory genes (XBP1,
MZB1, DERL3, PDIA4, HYOU1, RPN1, SYVN1, FKBP11, PRDM1) mean **0.1932** lognorm at 21.2%
detection; its 16 canonical B genes (MS4A1, CD79A/B, CD19, PAX5, EBF1, BANK1 …) mean
**0.0047** at 0.6%. A **41x ratio**. `XBP1` alone is **39%** of the whole core mean and is
detected in 64.1% of cells; `MS4A1` in **0.90%**. Acinar cells are the most secretory cells in
the body, so on this tissue the core reads the secretory program. **Split the core into its
two halves and score them separately — that settles the 26.5% in one run.**

**This is the k-instability.** Nonspatial k=20: group 7 wins `Exocrine` at z=3.597 with
`Immune_B` runner-up; group 2 wins `Immune_B` at z=3.539 with `Exocrine` runner-up. A z-gap of
**0.06** decides a 100%-vs-0% assignment, and at k=30 four groups flip. Two cores are reading
the same secretory signal and the ranking has no null of its own.

**Lineages are awarded on negative contrasts.** `nullgate.py:155` gates on `real > p95` with
no `real > 0` requirement. Nonspatial k=30 group 11 -> `Vascular` at **real = −0.0134**;
spatial k=30 group 12 -> `Myeloid` at **real = −0.0045**, 833 cells. Those groups read *lower*
on the awarded core than the rest of the tissue. Both are in `composition.csv` now.

**`Immune_T -> 0.0%` is a reporting defect, not a measurement.** `run_donor.composition`
(`run_donor.py:75-92`) builds rows only for lineages some group won, so a lineage winning
nothing has **no row**, and a reader turns that into a zero. Immune_T in fact clears its null
and is runner-up in six spatial groups (z up to 18.6). It was measured, then discarded by the
ranking. Emit a row for every lineage with an explicit status — `awarded` /
`cleared_but_outranked` / `no_group_cleared` / `no_gene_on_panel` — and `frac = NaN`, never 0.

**The pilot donor is also the wrong one.** 6539's PhenoCycler side has **124 peri-islet cells
over 541 islets** — its `Islet_N_exp20um` annotation was run on 3 of 539 islets, and it is the
only donor with no `Annotation (Tissue)` at all. 19 of 20 donors are >=98.8% ringed. So 6539
cannot anchor any peri-islet cross-modality comparison until QuPath is re-run on it.

## 3. Code built (untracked: `core/phenocycler/scripts/xenium_celltype/`)

| file | state |
|---|---|
| `nullgate.py` | **7 tests pass.** Annotates *groups*, awards a lineage only where the real core beats its size- and detection-matched null; groups clearing nothing → `Unavailable` with a reason |
| `xenio.py` | Dual-source loader (h5ad for 16 donors, bundle `cell_feature_matrix.h5` for the 4 brusko-only ones). Cross-validated on 6539: identical gene order, scores to 6e-7, **100% argmax agreement** over 20,000 cells |
| `run_kodama.R` | Both arms; records timings + KODAMA version to JSON |
| `run_donor.py` | PCA → KODAMA → KMeans → nullgate → composition with binomial CIs |
| `tests/test_nullgate.py` | 7 tests |

**Two design flaws its own tests caught** (both fixed, both pinned by a test):

1. The statistic had to become a **contrast** (in-group minus out-of-group). On a *level*, the
   brightest core wins every group while distinguishing none.
2. A **single group cannot be a contrast** — it now raises rather than falling back to the
   level, which would reintroduce the defect.

KODAMA **3.3** installed, pinned by commit. `smith6jt-cop/KODAMA` is an unmodified fork of
`tkcaccia/KODAMAlegacy` (0 ahead, 0 behind) — a pin, matching the DCT/UniFORM pattern here.

---

## 4. The pilot result — and why it is not yet trustworthy

Donor 6539, 50,000 cells, k=30, both arms. **The arms disagree enormously:**

| lineage | nonspatial | spatial | Δ |
|---|---|---|---|
| Exocrine | 32.0% | 51.9% | **19.9 pp** |
| Immune_B | 26.5% | 13.7% | **12.8 pp** |
| Unavailable | 15.4% | 6.1% | 9.3 pp |
| Myeloid | 2.5% | 5.9% | 3.4 pp |
| Immune_T | 2.7% | **0.0%** | 2.7 pp |
| Endocrine | 1.66% | 1.65% | 0.01 pp |

This vindicates running both arms: letting a cell's neighbourhood inform its own type call
changes the answer by up to 20 points, and **`Immune_T` vanishes entirely** in the spatial
arm. Only the `nonspatial` arm may carry a peri-islet spatial claim; using the spatial arm for
that would be circular.

**Do not report these numbers.** Two reasons, in order:

1. **The acinar no-axis problem (§2b) is unaddressed.** ~80% of the tissue has no core of its
   own, so it is being absorbed — most likely by `Exocrine` (the nearest, ductal) and possibly
   by `Immune_B`. Until that is handled, no share here means what its label says.
2. `Immune_B` at 26.5% is the flag to chase. **Not because B cells "should" be rare** — that
   reasoning is forbidden here and it is how the frequency rule gets broken — but because a
   58-gene core clearing its null on a quarter of the tissue is a mechanism worth
   understanding, and §2b supplies a candidate mechanism.

Endocrine at 1.66% is the one share stable across arms.

---

## 5. The cost model is UNRESOLVED, and I got it wrong three times

I quoted 25,000–62,000 core-hours. The maintainer rejected it as obviously miscalculated and
was right. **Three compounding errors, all in the same direction:**

1. **KODAMA parallelises over M, not over cells.** The M runs are the unit handed to
   `mclapply`. Every probe ran at **M=5 with `n.cores=16`**, so at most 5 cores could be busy
   — the 565% CPU I saw and misread as "poor parallel efficiency." It is not inefficiency;
   I asked for 5 units of work and offered 16 workers.
2. **The machine was contended.** Timings were taken while the peer session's `run_step1.py`
   (1719% CPU) and workflow agents' own probes (3930% CPU) were running. Relative results
   (PCs vs genes, back-to-back) survive; absolute seconds are upper bounds.
3. **I extrapolated an exponent from a contaminated three-point ladder.**

### What is actually measured (n=20,000, 50 PCs, idle-ish machine)

| operation | time |
|---|---|
| `Rnanoflann::nn` k=51 (landmarks=50) | 30.5 s |
| `Rnanoflann::nn` k=501 (landmarks=500) | 33.4 s |
| KODAMA M=1, ncomp=50, landmarks=500, 1 core | 77.8 s |
| KODAMA M=4, ncomp=50, landmarks=500, 4 cores | 219.6 s |
| KODAMA M=1, **ncomp=10** | 76.8 s |
| KODAMA M=1, **ncomp=50** | 76.8 s |

**Read these carefully — two hypotheses of mine died here:**

- **`ncomp` makes no difference** (76.8 s either way). I had predicted `ncomp=50` on a 50-PC
  input was degenerate and expensive. It is degenerate — but it is not the cost.
- **M-parallelism is not delivering.** M=4 on 4 cores took 219.6 s where M=1 on 1 core took
  77.8 s. Subtracting the ~33 s kNN, that is ~187 s for 4 parallel runs against ~45 s for one
  — i.e. **running essentially serially despite `n.cores=4`**. Cause unknown: `mclapply` fork
  copy cost, or memory-bandwidth bound. **This is the first thing to investigate**, because
  every cost projection depends on it.

The kNN is ~33 s and barely depends on k, so it is ~10% at n=20,000 and amortises over M
(built once, at source line 42; the parallel M-loop starts at line 217).

### Verified in the KODAMA 3.3 source (these are solid)

- `landmarks` **is** the principled subsample: `Xdata = data[landpoints,]`,
  `Tdata = data[-landpoints,]` — KODAMA optimises on landmarks and propagates to every
  remaining cell. **So `run_donor.py`'s `--max-cells 50000` default is wrong**: subsampling to
  50k and then setting `landmarks=500` subsamples twice.
- `neighbors = max(1L, floor(min(c(landmarks, nsample*0.75-1), 500)))` — the kNN graph is
  **hard-capped at 500**, so `landmarks` above 500 does not widen it. Upstream's
  `landmarks=10000` buys nothing there.
- `splitting = ifelse(nrow(data) < 40000, 100, 300)`. A design agent claimed my superlinear
  scaling was an artifact of this step. **I checked and it is not**: 5k→20k sits entirely
  within `splitting=100` and still gives exponent 1.41, and `splitting` drives
  `kmeans(Xdata, splitting)` on the *landmark* set (fixed size), so it cannot scale with n.

---

## 6. Next steps, in order

1. **Investigate why M-parallelism does not scale** (§5). Everything downstream is gated on it.
   Then re-time with `n.cores ≥ M` on an idle machine and build a cost model from
   measurement, not extrapolation.
2. **Solve the acinar no-axis problem** (§2b) before any composition is reported. Options: a
   data-driven cluster characterised by its own top genes and explicitly *not* given a
   taxonomy label; or `Unavailable, reason=no_axis_for_this_compartment`. Do not invent a core.
3. **Remove `--max-cells 50000` as the default** — it double-subsamples (§5).
4. **Decide per-donor vs joint embedding.** Joint makes clusters correspond across donors
   (which the donor-level comparison needs) but risks batch entering the embedding — and batch
   is the confound this study exists to break. Unresolved; the design workflow's answer has
   not been verified.
5. **HiPerGator.** Direct ssh from this workstation fails (`Permission denied
   (publickey,keyboard-interactive)`) — a key or the tunnel route is needed before anything can
   be submitted. Conventions to follow: `--account=maigan --qos=maigan-b`, work area
   `/blue/maigan/smith6jt`, R libs `/blue/maigan/smith6jt/R_libs`.
6. **Integration workstreams B/C/D** in the plan file remain untouched, including the
   contract/vocabulary rework and the crossed-batch de-confounding figure.

---

## 7. Traps

- **`S8 donor.py` cannot run at all today** — `run_donor` calls `require_current_pheno`, which
  raises unconditionally. The blocker is **PhenoCycler's** typing, not Xenium's. Do not name
  S8 as the consumer of this work until that is resolved.
- **The machine is shared** with the peer session, which is actively running heavy jobs on the
  same branch. Check `uptime` before timing anything.
- **This workstation is 32c/64t, 503 GB, RTX 4000 Ada 20 GB** — not the 16-core box my early
  framing implied. I passed `n.cores=16` by copying the benchmark's SLURM convention.
- **`/data/Xenium_Analysis/scripts/benchmark/`** already benchmarks KODAMA vs BANKSY/UTAG/MISTy
  — but for spatial **domains**, on different samples, with HiPerGator paths. It does not
  settle the cell-identity question and its "KODAMA is slowest" finding is explained by it
  feeding the gene matrix directly (PCs are 35× faster at identical settings).
- **20 post-Xenium re-stain qptiffs** are on HPG at
  `/orange/maigan/Xenium_T1D_PancPLN/Phenocycler/Batch<N>_<donor>_<serial>/Scan1/` — the
  directory names carry the donor ↔ serial mapping, which largely solves the
  `restain_sections.csv` problem flagged in the plan. Only 5 are local.

---

## Uncommitted state

```
 M core/phenocycler/config.ini
 M core/phenocycler/data/integration/xenium_paths.csv
 M core/phenocycler/docs/DEEPCELL_TYPES_REVIEW_2026-08-24.md
 M core/phenocycler/phenocycler/integration/manifest.py
 M core/phenocycler/tests/integration/test_integration_manifest.py
?? core/phenocycler/scripts/xenium_celltype/
```

`?? phenotype_workflow.csv` at the repo root is the **peer session's**, not mine.

Verification baseline: `pytest -q tests/` → **574 passed, 2 skipped**
(`geopandas`/`SimpleITK` absent; both are in the declared `integration` extra and should be
installed before the IoU and non-rigid registration paths are trusted).
`pytest -q scripts/xenium_celltype/tests` → **7 passed**.
