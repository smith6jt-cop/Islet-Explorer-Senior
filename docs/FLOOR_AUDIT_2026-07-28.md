# Raw input floor audit — 2026-07-28

> **Status: IN PROGRESS.** The audit instrument is being written. Findings below that are already
> established are marked as such; the per-marker recommendations are not written yet.

Entry point for the round that precedes this: [`GATE5_HANDOFF_2026-07-28.md`](GATE5_HANDOFF_2026-07-28.md).
This document supersedes that handoff's open questions 3 and 5, and corrects two of
[`GATE5_ROUND_2026-07-27.md`](GATE5_ROUND_2026-07-27.md)'s findings.

---

## Reversal: donor 6566 INS is NOT absent

**Maintainer, 2026-07-28: "6566 has clear INS expression."**

`data/restore_pair_validation/expression_reviews.csv` carried
`6566,INS,CONFIRMED_ABSENT`. That row is **removed**. It was wrong, and its own evidence field said it
could not yet be relied on:

> *"Provisional 2026-07-27, **from screen diagnostics NOT image review** … Recorded so the run does not
> halt on a biological absence. **MUST BE CONFIRMED FROM THE IMAGE** before any INS-dependent result is
> reported."*

**Maintainer, 2026-07-28 (second message): "There are definitely not over 4000 INS cells, maybe a couple
hundred at most for 6566."**

So **both** numbers on record are wrong, in opposite directions, and the truth is ~10²:

| quantity | value | verdict |
|---|---|---|
| `lineage_positive_n` shipped in `apply_manifest.csv` | **0** | wrong — INS is expressed |
| `threshold_positive_n` RESTORE called above the divisor | **4,135** | wrong — ~20× too high |
| maintainer, from the image | **a couple hundred** | ground truth |

The state_reason production recorded — *"target-arm fold 1.67 is below 2 and **image review confirms
absence**"* — contradicts the evidence field of the row that generated it. No image review had occurred.

**Both errors trace to the same collapsed floor.** Production's `linear_otsu` floor for 6566 INS is
**17.4**, which admits **412,459 cells — 23% of the whole tissue** — into the INS training arm. An arm
that is 23% of the tissue has nothing to separate, so the NNMF learns nothing (`class_fold` 3.4,
`separation` 0.071), the divisor that follows is meaningless, and it overcalled by ~20×. The technical
gate then caught the instability (`divisor_reproducibility` 0.690) and the expression review converted it
to a biological absence, swinging the answer to zero.

The audit's candidate methods on this donor-marker, which is now the round's one anchor with a known
answer:

| method | floor | cells called | separation | class fold |
|---|---|---|---|---|
| **linear_yen** | **745.3** | **473** | **0.645** | **125.0** |
| log_yen | 94.0 | 1,319 | 0.424 | 62.0 |
| linear_triangle | 79.3 | 1,542 | 0.375 | 54.0 |
| log_triangle | 55.1 | 4,730 | 0.162 | 21.1 |
| **linear_otsu** (shipped) | **17.4** | **412,459** | **0.071** | **3.4** |
| log_otsu | 7.4 | 1,207,764 | 0.038 | 4.3 |
| log_li | 5.2 | 1,425,222 | 0.027 | 5.9 |
| log_gmm_crossover | 4.6 | 1,471,823 | 0.024 | 6.6 |

The raw input floor is not itself the positivity call, so 473 is not directly comparable to "a couple
hundred". What matters is that **`separation` ranks the eight methods in exactly the order the
maintainer's ground truth implies**, on a donor-marker where the answer was known independently. That is
the frequency-free diagnostic working, and it is why the audit reports it.

**Three things follow.**

1. **The reasoning was a cohort-conformity argument, which this project's working rule forbids.** The
   evidence cited *"INS target_fold collapses 199/137 (AAB/ND) → 8.4 (T1D)"* — a comparison of this
   donor against other donors' values and against its clinical label, not against its own data.
2. **The mechanism was the input floor, not the biology.** 6566's INS floor is **17.4**, against
   191–2403 on the donors whose INS behaves, and **81.4% of its cells sit above it**. The arm separation
   collapsed because the floor collapsed.
3. **It removes the confound in the floor evidence.** A collapsed Otsu floor is *expected* where there is
   genuinely no bright mode, so beta-cell loss and floor failure could not be told apart from the floor
   alone. A donor with confirmed INS expression and a collapsed floor settles it: **the floor fails
   independently of whether the marker is expressed.**

This was the cohort's **only** biological-absence call in 480 donor-pairs, and it was selected by a
technical gate (`divisor_reproducibility = 0.690`) rather than by evidence of absence.

## The finding that prompted the audit

`double_high_fraction` — the statistic the 2026-07-28 exclusivity rule gates on — is defined as *"above
both markers' **raw input floors**"* (`restore_validation.py:2805-2810`): raw-MFI positivity, which the
project's working rule forbids as a definition of positivity. The same floors build the NNMF's training
arms (`balanced_arm_fit_indices:1354-1367`), so:

```
raw input floor → exclusive training arms → NNMF separator → Step-1 anchor → control set → Step-2 divisor → production call
```

A collapsed floor does not merely spoil a diagnostic. It spoils the divisor.

### Floor behaviour, measured from `screen_v12_production` (no new compute)

`linear_otsu` is the failing policy, and its failures are **donor-specific blowouts, not levels**:

| marker | floor method | % cells > floor (med / min / max) | floor spread | donors > 50% |
|---|---|---|---|---|
| **INS** | linear_otsu | 9.3 / 2.3 / **85.8** | **333×** | **7** |
| **SST** | linear_otsu | 7.7 / 3.6 / **90.3** | **38×** | **9** |
| **CD11b** | linear_otsu | **32.5** / 6.1 / 67.1 | **24.5×** | 2 |
| **B3TUBB** | linear_otsu | 3.7 / 2.1 / **55.9** | 5.7× | 1 |
| E-cadherin / Pan-CK / Ker8_18 / EpCAM | linear_otsu | 63–80% on every donor | 2.5–3.3× | 17–20 |
| all 12 `max(Otsu, Triangle)` markers | — | ≤ 21% | ≤ 8.9× | **0** |

The discriminator is each marker's **own donor-to-donor stability**, not an assumed abundance: the four
epithelial markers are high on every donor with tight spread, which is what a genuinely abundant marker
looks like. SST cannot be 3.6% positive on one donor and 90.3% on another.

### Corrections to `GATE5_ROUND_2026-07-27.md`

1. **The "SST/INS violate mutual exclusivity on 45%/35% of donors" finding is WITHDRAWN.** Corrected for
   the marginals the floors set, both are 20–100× *more* exclusive than chance:

   | pair | P(tgt⁺) | P(ref⁺) | expected | observed | obs/exp |
   |---|---|---|---|---|---|
   | `INS ← Pan_CK` | 0.093 | 0.907 | 0.0843 | 0.0004 | **0.01** |
   | `SST ← Pan_CK` | 0.077 | 0.924 | 0.0707 | 0.0022 | **0.05** |
   | `CD11b ← CD3e` | 0.959 | 0.043 | 0.0413 | 0.0015 | **0.04** |
   | `CD11b ← EpCAM` | 0.325 | 0.860 | 0.2744 | 0.1733 | **0.64** ← worst of 41 |

2. **Finding 3 (SST's bimodal `target_fold`) and Finding 7 (INS) share one mechanism with the above** —
   the floor collapse produces the bimodality, the low fold and the double-high blowout together. One
   cause, three symptoms.
3. **`CD11b ← EpCAM` survives as the worst pair**, but at 16× rather than 115×. CD11b's own floor spans
   24.5× across donors and calls a third of all cells positive, so whether CD11b's behaviour is its
   reference or its own floor **cannot be decided until the floor is fixed**.
4. **The exclusivity veto is withdrawn as a gate** (maintainer, 2026-07-28). Pair selection is either
   biological or algorithmic; ours is biological, and the manuscript's remedy for a biologically-exclusive
   pair that reads as correlated is explicit: *"if there exists non-specific staining in
   biologically-known mutually exclusive markers, the marker pair will be positively correlated. In this
   case, we can use the other mutually exclusive markers to infer background."* That is a staining-quality
   finding, not a pair-validity verdict. `CD11b ← EpCAM` stands as frozen on 2026-07-24; the
   double-positive columns remain as reported diagnostics.

### Other things this surfaced

- **Half the gate markers are thresholded on a different feature than the other half.**
  `qupath_compartment_full_20donor/` holds only **12** markers, so the 14 added in the Gate-4 expansion
  (INS, GCG, SST, CD34, SMA, Podoplanin, Keratin_5, CD79a, CD11c, Iba1, CD206, MPO, CD4, CD8) are
  thresholded on **whole-cell means alone**, while the original 12 use `Membrane+Cell`. A fallback, never
  a decision, never evaluated.
- **`divisor_statistic` defaults to `max` in code** (`config.py:172`) while the frozen `p99` lives only in
  the parent `config.ini:62`. A run that omits `--config config.ini` silently selects a policy that
  cannot complete.
- `control_fraction_le_1` is 0.9899–0.9900 on all 479 VALID rows — 1% by construction at p99, carrying no
  information. The informative quantity is `manuscript_divisor / canonical_divisor` (median 5.96×, max 200×).

---

# The audit

`scripts/gate5/audit_marker_floors.py` (new, committed — the 2026-07-21 instrument was ad-hoc and no
longer exists). **57 viable markers × 20 donors × 12 candidate floor methods = 13,644 rows**, 4.5 min.
Output in `data/floor_audit_v2/`: `floor_audit_per_donor.csv`, `floor_audit_summary.csv`,
`floor_classification.csv`, `floor_audit.pdf` (58 pages, one per marker) and `figures/<marker>.{png,csv}`.

It calls the **production** `raw_otsu_floor` / `raw_triangle_floor` rather than reimplementing them, on
the raw pre-REDSEA whole-cell mean over every canonical cell with a finite non-negative value — the
exact input `estimate_marker_input_floor` receives. Verified pair-invariant first: across the 500
(donor, marker) combinations in `screen_v12_production`, a marker's floor is bit-identical in every pair
it appears in.

**Both validation checks pass.**

- `linear_otsu` reproduces the floors production shipped on **280 donor-markers, max relative delta
  2.4e-16**. The audit and the screen are measuring the same thing.
- CD20 against `marker_floor_audit_v1`: `linear_otsu` 0.398 (v1: 0.494), `linear_triangle` 0.004
  (v1: 0.012). Same ordering, same order of magnitude. Exact agreement was never expected — v1 ran
  before `cell_qc` removed 6,303,013 objects that are not cells (20.04% of the population), and an
  Otsu threshold is a property of the population.

## The root cause: `raw_otsu_floor` clips at the 99.5th percentile

```python
upper = float(np.quantile(values, upper_quantile))   # upper_quantile = 0.995
clipped = np.clip(values, 0, upper)
histogram, edges = np.histogram(clipped, bins=bins, range=(0, upper))
```

The clip exists to stop one extreme value collapsing a 512-bin **linear** histogram into a single bin.
But **for a marker whose positive population is smaller than 0.5% of cells, the clip deletes the entire
positive population before Otsu sees it**, and the threshold necessarily lands inside the background.

| donor / marker | clipped floor (shipped) | → % positive | unclipped floor | → % positive |
|---|---|---|---|---|
| **6566 INS** | 17.4 | 23.2% | **864.2** | **0.02% — 403 cells** |
| 6380 INS | 13.3 | 22.5% | 1906.2 | 0.00% |
| 6493 SST | 50.6 | 25.8% | 756.0 | 0.20% |
| 6461 SST | 46.1 | 38.1% | 576.7 | 0.19% |
| 115 EpCAM | 616.5 | 50.5% | 604.6 | 53.1% |
| 6461 E-cadherin | 924.0 | 26.6% | 926.2 | 26.5% |

**403 cells on 6566 INS against the maintainer's independent image estimate of "a couple hundred at
most"** — versus 412,459 admitted by the shipped floor and 0 emitted by the withdrawn expression review.
That is the audit's one anchor with an externally known answer, and it lands on it.

The mechanism is specific: linear-space Otsu takes its leverage from the bright tail, which is exactly
what the clip removes. **Log-space binning is not the fix** — unclipped log-space Otsu puts 6566 INS
back at 1,215,101 cells, because in log space a rare bright population has no leverage at all.

## What it costs across the gate set

Classification in `floor_classification.csv`, by whether relaxing the clip changes the worst donor's
separation by more than 0.05:

**A. Clip-affected — 14 of 26** (`INS SST CD20 CD79a B3TUBB CD11c SMA GCG CD68 CD206 CD8 MPO Podoplanin
CD3e`). Worst-donor separation gains 0.06–0.41. INS 0.011 → 0.416, SST 0.073 → 0.434, B3TUBB
0.133 → 0.485.

**B. Clip-neutral — 9 of 26** (`E_cadherin EpCAM Pan_Cytokeratin Ker8_18 Vimentin CD31 CD34 CD11b
Iba1`). Gain ≤ 0.003. These are the abundant markers whose positive population is far above the clip.
**The four epithelial gate markers and CD31 are in this group, so the Gate-4 epithelial result is not
affected by this defect.**

**C. Unresolved at any tested floor — 3 of 26** (`CD4 Keratin_5 CD163`). No method reaches 0.30
worst-donor separation. These markers may simply not have a resolvable positive population in this
panel; that is a marker/staining question, not a method choice.

`max(Otsu, Triangle)` is **not** a blanket fix: it is the shipped policy for 8 of the 14 clip-affected
markers, and it is *worse* than plain unclipped Otsu on all of them, because taking the higher of two
thresholds cuts further into an already-truncated positive population.

## What the audit does NOT settle, and why there is no automatic recommendation

I built one and then deleted it, because neither candidate rule survives inspection:

- **Maximin `separation`** picks a floor calling **50.8%** of cells CD79a⁺, 37% CD20⁺, 37% CD4⁺. For a
  rare marker, splitting the background in half scores higher than the correct tail split, because
  between-class variance is maximised mid-distribution.
- **Requiring the floor to clear the background peak** then picks floors calling 0.01% for CD20, CD163
  and CD11b — near-empty.

So `separation` is a good failure *detector* and a bad *selector*. **The per-marker floor decision needs
the figures and your eye, which is what the plan promised.** The one thing the audit does settle on its
own is that the clip is a defect, and that removing it cannot hurt group B.

## The review surface

`data/floor_audit_v2/floor_audit.pdf`, worst consistency first. Each page, one marker:

- **left** — every donor's own log-intensity distribution as a ridgeline with **that donor's own
  production floor** drawn on it. Red where this donor's split is under half the strength of the same
  marker's typical donor. On the INS page the red floors sit visibly *inside* the grey background hump
  and the green ones far right on the beta-cell shoulder.
- **middle** — positive fraction by donor, one line per method, production emphasised.
- **right** — separation by donor and method, with the half-median line the red/green call is made against.

The headline statistic is **`separation_consistency` = min(separation) / median(separation) across
donors**: a marker measured against *itself*. It separates the two things an absolute cutoff confuses —
a floor that finds a real split on every donor even where the marker is abundant (EpCAM 0.92,
E-cadherin 0.89, Vimentin 0.90) from one that works on some donors and not others (INS 0.02, SST 0.16,
CD163 0.26). No frequency assumption anywhere.

## Also fixed this round

**`divisor_statistic` no longer defaults to a policy that cannot run.** It defaulted to `max`
(`config.py:172`, mirroring `PairValidationConfig` "so this indirection cannot change behaviour on its
own"), while the frozen production `p99` lives only in the **parent** repo's `config.ini`. A run
launched without `--config config.ini` therefore silently selected `max`, which halts `MODEL_UNSTABLE`
on 115 `Vimentin←E_cadherin`, 6436 `B3TUBB←EpCAM` and 6591 `CD68←E_cadherin`. Both policy keys are now
unset by default and `config_from_pipeline` refuses to guess, naming the file to fix. The test that
asserted the mirroring is replaced by one asserting the refusal, plus a test that the parent
`config.ini` really does carry the frozen values. **254 tests pass** (was 252).

**The three orphaned Gate-5 tables are scripted and reproduced.**
`scripts/gate5/gate_counterfactuals.py` (new) regenerates `counterfactuals.csv`, `followup.csv` and
`mes_loss.csv`, which backed Findings 2, 4 and 8 but had no producing script anywhere in the repo. It
calls `lineage.assign_donor` with `COMPARTMENT_GATES` temporarily swapped — the technique
`scripts/gold_standard/score.py` already uses — so a counterfactual arm cannot drift from what
production would do. Every published headline reproduces to the decimal over all 25,151,028 cells:
CD11b⁺ epithelial-marker⁺ **81.1%**, and outside islets **80.8%**, no other immune marker **78.6%**;
of CD31⁺ also CD34⁺ **85.8%**; of CD34⁺ also CD31⁺ **39.7%**; CD34⁺CD31⁻ **10.5%** of all cells,
**44.5%** Vimentin⁺, **1.0%** in islet core. The Finding-4 counterfactual table reproduces exactly.

Reproducing it recovered an undocumented definition: **"epithelial-marker positive" in those tables
means the 5 STRUCTURAL markers** — the shipped `COMPARTMENT_GATES["Epithelial"]` **minus the three
hormones**. Verified by matching donor 6591 to the cell (`cd11b_epi` 13,454; `cd34_not_cd31_epi`
26,123). It is a defensible choice for these questions — "is this CD11b⁺ cell also cytokeratin⁺" is a
different claim from "is it INS⁺" — but it was never recorded. The script now emits **both** definitions.
The difference is 0.3% for CD11b, so the 81.1% headline stands, and 5% for `cd34_not_cd31_epi`.

## What this means for the open decisions

- **CD11b is clip-neutral** (gain 0.003, consistency 0.63, no donor collapsed). Its behaviour is *not*
  a floor artefact, so the `CD11b ← EpCAM` question stands where it was — and since the exclusivity
  veto is withdrawn, it stands as you froze it on 2026-07-24.
- **CD34 is clip-neutral** (gain 0.000, consistency 0.90) and among the cleanest markers in the panel.
  Removing it from the Endothelial gate remains purely a lineage decision on the evidence already
  measured, unaffected by anything here.
- **The Gate-4 epithelial result is unaffected** — all four epithelial gate markers are group B.
- **The 14 group-A markers cannot be trusted at their current floors on their collapsed donors**, and
  they include 8 of the 10 Immune gate markers. Any Immune claim rests on floors this audit flags.

## Next, and not started

Changing a floor method re-derives the NNMF training arms and therefore invalidates the screen, the
reference dossier, `restore_apply`, `lineage` and the Gate-5 dashboard. **Nothing downstream has been
touched.** The decision is per marker, from `floor_audit.pdf`.
