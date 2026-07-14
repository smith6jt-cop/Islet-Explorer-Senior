# HANDOFF — partner-low thresholds, proliferation binarization, PhenotypingPlan reconciliation

_Written 2026-07-14. Branch (both parent + submodule): **`redsea-restore-curated-5compartment`**. Nothing new pushed._

## TL;DR
The REDSEA/RESTORE **migration is committed** (submodule `123d916`, parent `a667ee0`). On top of that, this
session added **new, UNCOMMITTED work** in the `Phenocycler_Analysis` submodule working tree (375 insertions /
7 files) + parent `config.ini`: a **partner-low endocrine threshold method**, **proliferation (Ki67/PCNA)
binarization**, a reusable **`plot_pair_cloud` figure**, and the **`PhenotypingPlan.txt` reconciliation**
(Collagen→mask). **60 tests pass.** Everything is **default-off / additive** — the committed pipeline behavior
is unchanged until a knob is set. Several **science decisions are open** (below) before this gets committed.

## Git state
- **Submodule `Phenocycler_Analysis/`** — branch `redsea-restore-curated-5compartment`, HEAD `123d916`.
  **Uncommitted** (`git diff HEAD`): `restore.py`(+174), `spatial.py`(+107), `tests/test_diagnostics.py`(+67),
  `marker_taxonomy.py`(+19), `config.py`(+15), `diagnostics.py`(+13), `pipeline.py`(+7). No untracked files.
- **Parent repo** — branch `redsea-restore-curated-5compartment`, HEAD `a667ee0`.
  Uncommitted: `Phenocycler_Analysis` (submodule dirty pointer), `config.ini` (added `partner_low_q`,
  `proliferation_k`), untracked `PhenotypingPlan.txt` (the user's spec — keep) + `figures/*.png` (below).
- Snapshot of the migration's pre-commit state is still in the scratchpad (`phenocycler_worktree.tgz`).

---

## What's UNCOMMITTED this session (ready to commit, pending the open decisions)

### 1. Partner-low endocrine threshold method (`restore.py`, `config.py`, `config.ini`, `diagnostics.py`)
Idea (from the user): threshold a mutually-exclusive marker using cells that are **LOW for its partner**.
- `bimodal_thresh(v, ...)` — 2-GMM crossover on `log10(nonzero+1)`, `+inf` if not bimodal (shared core).
- `partner_low_thresh(T, R, partner_low_q, islet=..., ...)` — conditions on reference-LOW (+ islet) cells,
  delegates to `bimodal_thresh`. **Endocrine-only** (INS/GCG/SST): wired in `patch_restore_idx_floor.wrapper`
  guarded by `xlabel in ENDOCRINE_MARKERS`; non-endocrine fall through to the normal `neg_stat`.
- `cell_region` is threaded through `build_restore_input(..., region_dir=cfg.cells_dir)` **only when
  `partner_low_q>0`**, so the wrapper can restrict to islet (`core`+`peri`) cells.
- Config knob **`restore_partner_low_q` (default 0.0 = off)**; `--partner-low-q` CLI.
- `r_otsu` was **relocated** from `diagnostics.py` into `restore.py` (re-exported) to avoid a circular import.

**Result (verified on 6374 ND + 6380 β-loss):**
- Bidirectional works — `MARKER_PAIRS` has both `[INS,GCG]` and `[GCG,INS]`, so partner-low computes **INS via
  GCG-low (=398)** AND **GCG via INS-low (=491)** on 6374. (SST is one-way, `[SST,INS]` only.)
- 6380 (β-loss) → **ABSENT** (correct — GCG-low INS is unimodal, no β mode).
- Extra-islet endocrine ARE identified: the islet cut is **applied to all cells**, so 6374 gets **3,405
  extra-islet INS⁺ + 2,164 extra-islet GCG⁺** (single cells + a clear α cluster — see figure).

### 2. Proliferation (Ki67/PCNA) binarization — GRADED, not bimodal
- `logmean_ksigma(v, k=3.0)` — robust per-donor upper-tail cutoff = `mean + k*std` of `log10(nonzero+1)` →
  raw MFI. (Ki67/PCNA have a zero spike + one broad mode + a tail — **no valley**, so a bimodal/Otsu cut
  gives a nonsense ~70% positive.) **User chose log mean+3σ.**
- `run_restore_proliferation(cfg, k=...)` applies it per (donor, marker) → writes `{m}_pos/_norm/_log2r` to
  `data/restore_gated_proliferation/` + `restore_thresholds_proliferation.csv`.
- `--proliferation` CLI; a **`proliferation` pipeline stage** (added to `STAGES`/`ORDER` after `restore`).
- `PROLIFERATION = ["Ki67","PCNA"]` in `marker_taxonomy.py`; config knob **`restore_proliferation_k` (=3.0)**.
- **Outcome:** Ki67 sensible (6374 0.10%⁺ thr 348; 6380 1.38%⁺ thr 579). **PCNA is problematic at k=3** —
  6374 PCNA → ~0 (1 cell, thr 3039), because PCNA is broadly expressed (high G1 baseline) so mean+3σ overshoots.

### 3. `plot_pair_cloud` (D1 figure) — `spatial.py`
Reusable 6-panel INS↔GCG figure: before/after REDSEA + Pearson r, **colored by parent object (core/peri/
tissue)**, current RESTORE thresholds (quadrants), islet-only, and the **partner-low panel** (islet reference-
low target histogram + GMM crossover vs current cut vs Otsu valley). Uses **DuckDB** to join cells+cells_redsea
(a pandas string-key merge on 2.5M rows hangs >4 min; DuckDB does it in 0.8s).

### 4. `PhenotypingPlan.txt` reconciliation (`marker_taxonomy.py`)
The user's `PhenotypingPlan.txt` (repo root) matches the code **except**: **IAPP** (plan includes `[IAPP,GCG]`;
code keeps it in `EXCLUDED` as a failed marker) — user confirmed **keep IAPP + DAPI OUT**. Done: **Collagen_IV
→ new `MASK` list**, removed from `TYPE` (out of the identity heatmap, used as a spatial mask); `heatmap_markers`
subtracts `MASK`. `FUNCTIONAL_PAIRS` (12) all match. Leave-out 6 / not-lineage-clean 8 all handled.

---

## OPEN DECISIONS (resolve these next)
1. **Partner-low activation & direction.** On 6374, partner-low INS = **398 (specific, 32% islet β)** vs the
   current `mean3sd` = **29 (sensitive, 82%)**. Which is biologically right? The user's earlier "~0.25–0.5 on
   both axes" is far below both — likely a different scale; the figures show the real valley at log10 ~2.4–2.6.
   Default stays **off** (`restore_partner_low_q=0.0`) until decided.
2. **PCNA handling** — lower `proliferation_k` for PCNA (per-marker k?), drop PCNA from `PROLIFERATION` and keep
   Ki67 only, or accept ~0. (Ki67 at k=3 is good.)
3. **SST bidirectional** — add `[INS,SST]`/`[SST,GCG]`/`[GCG,SST]` reciprocals to `MARKER_PAIRS` if wanted.
4. **D2 — REDSEA double-positive diagnostic (NOT STARTED).** Classify the ~17% INS⁺GCG⁺ "double-positive"
   cells (β↔α boundary/segmentation artifact vs genuine) via the contact matrix + ROI renders; bounded
   `--alpha-per-channel`/gate sweep with an over-correction guard. Diagnose + recommend only (user's choice);
   do NOT change production REDSEA. Design is in the (superseded) plan file — see below.
5. **Commit the uncommitted work** once the above are settled; then bump the parent submodule pointer.
6. **Full 22-donor recompute** is deferred (all validation so far is on 6374 + 6380).

## Figures (all untracked in `figures/`)
- `6374_INS_GCG_cloud.png`, `6380_INS_GCG_cloud.png` — the 6-panel pair cloud (region + partner-low).
- `6374_extra_islet_endocrine.png` — extra-islet β/α spatial (single cells + clusters).
- `proliferation_calibration.png` — Ki67/PCNA distributions + candidate cutoff rules.
- `6374_ins_gcg_cloud_thresholds.png` — the original 4-panel that started this.
- `eval_item1_negstat.png`, `eval_item2_presence_6380.png` — the 3-knob evaluation (all knobs stay off).

## How to resume / run
```bash
conda activate scvi-env            # env has phenocycler (editable), spams, duckdb
cd /home/smith6jt/IO60panc2nd/Islet-Explorer-Senior
CONFIG=$PWD/config.ini             # ALWAYS pass --config (default loads the submodule template)
# tests
cd Phenocycler_Analysis && python -m pytest tests/ -q     # 60 pass
# proliferation binarization (log mean+k*sigma)
python -m phenocycler.restore --config $CONFIG --proliferation --donors 6374 6380
# partner-low eval (endocrine, islet, opt-in) — thresholds only
python -m phenocycler.restore --config $CONFIG --donors 6374 6380 --partner-low-q 0.5 --skip-apply
# regenerate the pair-cloud figure
python -c "import matplotlib;matplotlib.use('Agg');from phenocycler.config import load_config;\
from phenocycler.spatial import plot_pair_cloud;cfg=load_config('config.ini');\
plot_pair_cloud(cfg,'6374','INS','GCG',thr_target=29.4,thr_ref=70.4)"
```
Per-stage CLIs take `--donors`; the top-level `python -m phenocycler.pipeline` has **no donor filter** (would
recompute all 22). `data/cells/` + `data/cells_redsea/` exist for all 22 donors.

## Gotchas learned this session
- **`df.T` is the DataFrame transpose**, not the column `"T"` — on 2.5M rows it silently hangs (object-dtype).
  Use `df["T"]`. (This caused a "figure never renders" bug.)
- **Join millions of cells with DuckDB, not a pandas string-key `.merge`/`.join`** (0.8s vs >4 min hang).
- **`conda run` + a pipe buffers stdout** — `print(flush=True)` still won't appear until exit; write timing to
  a file (not `| grep`) when debugging, and use `--no-capture-output`.
- `partner_low_thresh`/`bimodal_thresh` return **`+inf`** for absent/unimodal → `_apply_one_donor` skips that
  (image, marker) cleanly (no `_pos` column) — intended (β-loss INS, unimodal proliferation).
- The old plan file `~/.claude/plans/the-reason-for-locating-tranquil-cosmos.md` describes D1/D2/D3 (D2 design
  still useful); it predates the partner-low/proliferation findings here.
