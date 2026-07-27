# Anatomy check — the Epithelial gate, not (only) the divisor

Round report, 2026-07-27. Reads on from
[`HANDOFF_2026-07-27_qc_and_v10.md`](HANDOFF_2026-07-27_qc_and_v10.md). Figures in
`figures/anatomy_check/`; every figure writes its own source CSV beside the PNG.

> `figures/*` is gitignored, so the PNGs and CSVs are **local-only** — regenerate them with
> `python scripts/anatomy_check/figures.py`. Every number this report relies on is in its tables, so
> the conclusions do not depend on the files being present.

**Nothing here is canonical.** Every arm ran in an isolated data root under `data/_arm/`; canonical
`data/` was not written. No divisor is frozen, no gate is changed in production, Gate 4 is not re-run.

---

## The question this round asked

The handoff left one open decision — which divisor policy — and named the anatomy check as the arbiter
rather than another threshold sweep. Before running it, three measurements changed what the check
needed to test:

1. **Production was wired to the worst arm.** `restore_apply.run_apply` accepted a
   `PairValidationConfig` but nothing supplied one, so a production apply always ran the v10 dataclass
   defaults `reference_only` + `max` — the arm that calls E-cadherin on 4.3–18.9% of cells, *worse*
   than the v9 arm the Gate-4 pilot was rejected on. Fixed (Phase 1); the knobs are now `config.ini`
   keys with CLI overrides, and the defaults still reproduce the previous behaviour exactly.

2. **One global divisor statistic moves every marker at once.** Median call rate over the 7 pilot
   donors, `ref_only/max` → `ref_only/p99`: E-cadherin 7.7 → 48.3%, but CD3e 0.08 → 5.5%, CD68
   0.14 → 5.1%, Vimentin 0.47 → 18.6%, B3TUBB 0.02 → 2.3%. Tuning the statistic for E-cadherin
   inflates every sparse marker by the same operation (figure D).

3. **No single E-cadherin threshold can cover both epithelial populations this tissue has.** Cohort
   medians over all 20 eligible donors, REDSEA values joined to the QuPath `cell_region`:

   | region | n | E-cad | EpCAM | Ker8/18 | Pan-CK | INS | GCG |
   |---|---|---|---|---|---|---|---|
   | islet core | 656,633 | **630** | **656** | 326 | 175 | 407 | 1004 |
   | peri-islet | 680,581 | 1004 | 541 | 1049 | 454 | 23 | 138 |
   | tissue | 23,813,814 | **762** | **516** | 887 | 367 | 10 | 71 |

   Endocrine sits *below* acinar on E-cadherin and *above* it on EpCAM — opposite polarity across the
   one boundary that matters. INS 407 / GCG 1004 in core against 10 / 71 in tissue confirms the
   `Islet_N` annotations are trustworthy anatomy, independent of any marker threshold.

So the round tested two hypotheses together, since they trade against each other: **H1** the divisor is
too strict, **H2** a single-marker Epithelial gate cannot work here.

---

## Method

Five divisor arms × two Epithelial gate designs × the 7 Gate-3 pilot donors (115, 6380, 6436, 6450,
6476, 6539, 6591).

- **Gate designs.** `single` = production today, `Epithelial = E_cadherin`. `union` =
  `E_cadherin ∪ Pan_Cytokeratin ∪ Ker8_18 ∪ EpCAM`. `lineage.py:146-161` already ORs across a gate's
  marker list — exactly as Immune does over CD3e/CD68/CD11b — so this is a config change, not new logic.
- **The three extra markers were run through the same frozen method** as production targets
  (`PAIR_SETS["epithelial-comparator"]`, each against CD31 and Vimentin), reusing
  `evaluate_locked_pair` and the tested per-cell mapper `restore_apply.cell_calls_from_evaluation`.
  They are **diagnostic only** — not in `ACCEPTED_PAIRS`, no production call.
- **Two controls on the diagnostic path itself.** The comparator's `E_cadherin ← CD31` divisor is
  **bit-identical to the production apply's** on every arm and donor (max |Δ| = 0, n = 7 per arm), and
  every comparator pair returned `VALID` on every arm and donor (14 rows per target per arm) — so the
  union gate inherits no technical failure and the comparator is not doing anything different from
  production.

Reproduce: `python scripts/anatomy_check/run_arms.py` then `python scripts/anatomy_check/figures.py`.

---

## Result 1 — one arm is eliminated on technical grounds

`reference_only/max` **cannot produce a production run at all.** `restore_apply` halted on donor 115
`Vimentin ← E_cadherin` with `MODEL_UNSTABLE`; the 5% screen had already flagged 3 of its 49
donor-pairs the same way. This is the production guard behaving correctly, and it removes the arm from
consideration before any biology is argued. It is also the arm production has been silently configured
to use since Gate 4 shipped.

## Result 2 — the gate is at least as large a lever as the divisor

**Islet-core cells** — endocrine by QuPath annotation, hormone-validated, unambiguously epithelial —
min–max across the 7 donors (figure B; per-donor points, not just the median):

| arm | Epithelial, single gate | Epithelial, union gate | `Other`, single | `Other`, union | `Immune` |
|---|---|---|---|---|---|
| v9 (`reference_and_target`/max) | **1–6%** | 16–61% | **82–94%** | 24–79% | 0–9% |
| `reference_only`/p99.9 | 2–10% | 7–55% | 76–90% | 32–85% | 1–4% |
| `reference_only`/p99 | 11–44% | **28–72%** | 31–76% | **6–60%** | 3–11% |
| `reference_only`/p95 | 35–71% | 56–84% | 2–47% | 0–29% | **7–31%** |

Read the first row: under the gate and divisor **production is currently wired to**, 82–94% of
indisputably epithelial islet cells are classified `Other`. Holding the divisor fixed and changing only
the gate moves islet-core Epithelial from 1–6% to 16–61%; holding the gate fixed and changing only the
divisor to p99 moves it from 1–6% to 11–44%. **The two levers are comparable in size and they compose.**

**H2 is supported. The single-marker Epithelial gate is a real defect, independent of the divisor.**

The last row is the compensating-error trade made explicit, and it is why the epithelial evidence
cannot settle the divisor on its own: p95 gives the best epithelial purity of any arm *and* puts
7–31% of islet-core cells in `Immune`. Whether that is insulitis or over-calling is not answerable
from any internal check — it needs the gold standard.

## Result 3 — but the union is not a complete answer either

Figure G, cohort means over the 7 donors:

| arm | union adds (% of all cells) | of which ≥2 markers agree | 0 of 4 positive | 4 of 4 positive |
|---|---|---|---|---|
| v9 | 12.7 | 3.8 (30%) | 65.9 | 5.0 |
| p99.9 | 7.1 | 2.9 (41%) | 64.7 | 3.8 |
| p99 | 11.3 | 4.0 (35%) | 40.4 | 14.3 |
| p95 | 11.6 | 6.9 (59%) | 22.8 | 34.9 |

So roughly **two thirds of what the union adds rests on a single marker**. The union is a genuine
recovery, not a licence to treat any single-marker positive as epithelium, and at the best viable arm
it still leaves 40% of cells negative for all four epithelial markers.

The pairwise Jaccard (figure G, right; mean of 7 donors at `reference_only`/p99) explains *why* it
works and names the marker doing the work: the three keratin-family markers agree closely with each
other (E-cad↔Ker8/18 0.72, Pan-CK↔Ker8/18 0.75, E-cad↔Pan-CK 0.64) while **EpCAM is the outlier at
0.31–0.37 against all three** — it marks a substantially different cell set. That is exactly the
polarity seen in the cohort medians above (EpCAM higher in islet core, the keratins lower), i.e. EpCAM
is recovering endocrine epithelium the keratins miss. The union is not averaging noise; it is adding a
biologically distinct epithelial population.

Figure F, per-marker distributions on donor 6539, shows the same thing directly and independently of
any threshold: islet-core EpCAM is shifted **right** of tissue, while islet-core Ker8/18 and
Pan-Cytokeratin are sharply peaked near zero. Endocrine is EpCAM-high / keratin-low; acinar is the
reverse. One marker cannot span that.

**A third, unrelated reason the union is more robust.** REDSEA clamps negative results to zero, and it
does so unevenly across these four markers. Cohort-wide over the 20 eligible donors, the fraction of
cells whose REDSEA value is exactly 0: E-cadherin **8.0%** of tissue cells (3.1% of islet core),
Pan-CK 4.8%, Ker8/18 4.0%, EpCAM 0.5% — but **only 0.35% are zero on all four**. A single-marker gate
inherits its marker's clamping rate as an irreducible floor of false negatives; a union of four does
not.

## Result 4 — a positive biological control the pipeline was not tuned to produce

Figure C. Endothelium is **enriched inside islets** (8–13% of islet-core cells vs 2–6% in tissue) and
dips at 5–10 µm before recovering. Islets are densely vascularised, so this is the correct answer, and
no threshold in this pipeline was tuned toward it. It is the strongest evidence so far that the
endothelial call is sound — and the contrast is what makes the epithelial result damning: on the same
axis, in the same run, Epithelium is at its *lowest* inside islets.

## Result 5 — the union raises spatial coherence; a stronger claim I first drew from figure E does not survive

Figure E, lower row: in a 350 µm window of solid acinar tissue, `Epithelial` and `Other` are finely
interleaved cell by cell. My first reading was that this is the signature of a threshold cutting
through one homogeneous population. **Figure H, which quantifies exactly that, does not support it.**

For each cell, take the share of its 12 nearest neighbours in the same compartment and rescale it as
`(observed − marginal) / (1 − marginal)`, so 0 is random labelling and 1 is perfect segregation. The
null comes from that donor's own composition, so nothing is assumed.

| compartment | v9 single → union | p99 single → union |
|---|---|---|
| Epithelial | 0.130 → 0.146 | 0.140 → 0.179 |
| Immune | 0.074 → 0.074 | 0.096 → 0.096 |
| Endothelial | 0.046 → 0.047 | 0.059 → 0.061 |
| `Other` | 0.073 → 0.070 | 0.045 → 0.045 |

Every compartment scores 0.02–0.20. Single-cell neighbourhoods in pancreas are genuinely mixed —
capillaries, immune cells and stroma thread through every acinus — so low absolute values are expected
and are not evidence of anything. Within that, **Epithelial is among the *more* structured**, not the
less, and **the union raises it at every arm**. `Other` is the least structured everywhere, which is
what a residual bucket of heterogeneous leftovers should look like.

Two things follow. The union's added cells are not scattered noise — they increase spatial coherence,
independently corroborating Result 3. And figure E should be read as illustration, not as evidence
that the epithelial call is a threshold artifact.

> Method note worth carrying forward: the first version of figure H used the ratio `observed/marginal`,
> whose ceiling is `1/marginal`. On that scale Neural (0.5% of cells) scored 26 against Epithelial's
> 1.5 and the conclusion inverted. Any neighbour-concordance statistic compared across compartments of
> different size has to be rescaled to a common range first.

---

## What this does and does not settle

**Settled.** `reference_only`/max is not a candidate. The single-marker Epithelial gate is a defect on
its own terms, measured against anatomy rather than against an expected frequency, and the union gate
both recovers islet-core epithelium (Result 2) and raises its spatial coherence (Result 5). The Gate-4
wiring reproduces the frozen method exactly (comparator ≡ production divisor, max |Δ| = 0).

**Withdrawn.** That the epithelial call is spatially structureless — figure H measures it and it is
not. See Result 5.

**Not settled.** The divisor policy. The union at `reference_only`/p99 gives the best islet-core purity
of anything tested, but p99 also multiplies every sparse immune marker 5–60× (figure D), and *that*
cannot be adjudicated by the epithelial evidence here. The per-donor spread stays wide (islet-core
Epithelial 28–72% under the best configuration), which is itself a finding: some donors are much worse
than others and the cohort median hides it.

**Still unmeasured.** Accuracy. Everything above is internal consistency plus anatomy; none of it says
what fraction of cells called `Immune` are immune cells, which is the question the study rests on. The
gold-standard bundle built this round (`data/gold_standard_v1/`, protocol in
`scripts/gold_standard/README.md`) is the instrument for that and needs annotator time, not code.

---

## What changed in code (Phase 1 — required regardless of the outcome)

| change | why |
|---|---|
| `[restore] control_definition` / `divisor_statistic` in `config.ini` → `PipelineConfig` → `restore_apply.run_apply`, plus CLI flags on `restore_apply.main()` | the v10 knobs were reachable only from `restore_validation evaluate-locked`, so the arm under review and the arm production ran could not be made to agree. Defaults still reproduce prior behaviour exactly; a bad value now fails before any donor is evaluated. |
| `EPITHELIAL_TARGET_COMPARATOR_PAIRS` + `PAIR_SETS["epithelial-comparator"]` | Pan-CK / Ker8/18 / EpCAM as **targets** (they were references only), each against CD31 and Vimentin. Adds no new QuPath marker, so the existing compartment sample answers it with no re-export. `ACCEPTED_PAIRS` untouched. |
| `cell_qc` added as pipeline stage 0 | a documented pre-estimation step that existed only as `python -m phenocycler.cell_qc`, so a clean re-run silently skipped it. Guarded to refuse a destructive re-run over already-filtered cells. |
| `Phenocycler_Analysis/scripts/run_divisor_sweep.sh` | the `sweep_qc` arms were ad-hoc CLI invocations recorded nowhere; `grep -rIn sweep_qc` over the repo returned nothing. Used it to finish the aborted p95 arm. |
| `scripts/build_acquisition_batch_table.py` → `data/donor_acquisition_batch.csv` | acquisition batch is a Gate-5 covariate that existed in no schema. Reproduces Fisher p = 0.0147 and adds the panel-variant split the audit did not record: 13 donors CD38, 4 b-Catenin1, **3 with neither (58 channels: 6523, 6566, 6591)**. |
| 227 tests (was 217) | new coverage for the config threading, the policy override precedence, the invalid-value guard, and the comparator pair-set invariants. |

Also corrected, each verified against the artifact rather than against another document: the
`cohort_manifest.json` cell count (31,454,041 → 25,151,028), the rebuild plan's "Gate 2 of 7" banner
and its Gate-6 cell count, `README.md`'s KINTSUGI claim, `restore_validation.py`'s "not yet wired into
lineage.py" comment, the `pipeline.py` stage-order docstring, `scripts/CLAUDE.md`'s r≈0.998 and
23.3M-cell claims, and `data/restore_pair_validation/README.txt`. The **311,340 vs 141,938 duplicate
count** is reconciled in the handoff: 311,340 cells in duplicate groups = 155,670 physical duplicates,
of which 141,938 were attributed to `duplicate_detection` and 13,732 to an earlier predicate that
caught the loser first.

## Recommended next step

Promote Pan-Cytokeratin, Ker8/18 and EpCAM to production RESTORE targets through their own Gate-1→3
validation — reference selection under the pure-proportion rule (both CD31 and Vimentin are already
screened here), then a full-cell stress pilot — and decide the divisor **after** that, because the
divisor question changes shape once Epithelial is a union. Do not freeze a divisor to compensate for a
gate that is about to change.

---

# Follow-on — Gate-2 reference selection, re-run on QC-clean data

The v1 reference dossier was stale: `reference_to_target_ratio` is computed from the negative-control
population, which the `control_definition` fix redefines, and the cell set changed underneath it. It
has been regenerated.

**Screen:** `data/restore_pair_validation/screen_v11_qc/ctrl-reference_only_div-p99/` — 27 pairs × 20
eligible donors = 540 rows, full cells, post-QC, `control_definition = reference_only`. 140 VALID, 398
`PAIR_REVIEW_PENDING` (the un-reviewed candidates, as expected), 2 `MODEL_UNSTABLE`. The decision
columns (`reference_to_target_ratio`, `candidate_full_maximum`) do not depend on `divisor_statistic`,
so one screen settles the reference question for any divisor policy.
**Dossier:** `data/restore_pair_validation/reference_selection_v2/`.

## Result A — every frozen reference survives

| pair | v1 ratio | v2 ratio | v1 → v2 |
|---|---|---|---|
| B3TUBB ← EpCAM | 64.1 | 52.9 | accepted → accepted |
| CD3e ← E_cadherin | 22.2 | 28.8 | accepted → accepted |
| CD68 ← E_cadherin | 7.4 | 5.9 | accepted → accepted |
| CD31 ← E_cadherin | 3.64 | 3.67 | accepted → accepted |
| CD11b ← EpCAM | 3.66 | 1.92 | accepted → accepted |
| Vimentin ← E_cadherin | 1.38 | **1.04** | accepted → accepted (marginal) |
| E_cadherin ← CD31 | 0.19 | **0.08** | accepted by exception → same |
| **CD11b ← E_cadherin** | 1.50 | **1.00** | comparator → **rejected** |
| **CD11b ← Ker8_18** | 1.12 | **0.90** | comparator → **rejected** |

Only two rows change disposition, and neither is a production choice. Two things are worth recording:

- **CD11b's original baseline reference now fails.** `CD11b ← E_cadherin` falls to 1.00 and is
  undersized on 50% of donors. The maintainer's 2026-07-24 override to EpCAM was made on other grounds
  (the baseline was "marginal"); under the corrected control definition it is the **only** passing
  reference CD11b has. The override was not merely defensible, it was necessary.
- **`Vimentin ← E_cadherin` is now the marginal one** — 1.04, undersized on 40% of donors. It is the
  next pair that would flip if anything upstream moves again. Worth watching, not acting on.

## Result B — the proportion rule has no admissible answer for the epithelial targets

All three promotion candidates come back `open`, and so does the already-accepted E-cadherin:

| target | best reference | median ratio | undersized on |
|---|---|---|---|
| Pan_Cytokeratin | Vimentin | 0.18 | 100% of donors |
| Ker8_18 | Vimentin | 0.18 | 100% of donors |
| E_cadherin | Vimentin | 0.16 | 100% of donors |
| EpCAM | Vimentin | 0.08 | 100% of donors |

Every candidate reference is undersized on every donor. **This is a property of the targets, not of
the references.** A marker positive on 25–50% of cells cannot have a *larger* negative control unless
some other marker is positive on more, and in pancreas none is. The rule was designed for the immune
targets, where it does real work — it correctly rejects `CD11b ← CD3e` (0.04) and `CD68 ← CD3e` (0.19)
— and the maintainer already made this exception once, for `E_cadherin ← CD31`. The screen now shows
that exception is not a special case: **it is the general situation for every majority-positive
structural target.**

What still discriminates for these pairs: `reference_fold_median` is 19–37× for all eight epithelial
rows, i.e. the arms are strongly separated and the fits are sound. The absolute control sizes are
large. Only the *relative* size criterion fails.

**This is a maintainer decision, not a computation.** Either the proportion rule gets a stated scope
(sparse targets only, with arm separation and absolute control size governing structural targets), or
each epithelial target gets an individually recorded exception like E-cadherin's.

## Two logic fixes this exposed

- `divisor_vs_consensus` was computed against the consensus of *rule-passing* references. With none
  passing, that set is empty and the column was `NaN` for exactly the targets whose reference choice
  most needs an independent cross-check. It now falls back to all of the target's references and
  reports `consensus_over_all_references` so the reader knows which basis was used.
- The `open` rationale said only "no candidate reference satisfies reference≥target frequency", which
  sends a reviewer looking for a better reference that cannot exist. It now names the structural cause
  when every candidate is undersized on every donor — keyed on that measured pattern, not on a cutoff.

233 tests.
