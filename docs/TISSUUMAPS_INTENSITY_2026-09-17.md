# TissUUmaps pixel intensity: what was wrong, and what it is now

2026-09-17. Scope: the reported defect "TissUUmaps shows the wrong pixel intensity distributions".
Everything below was measured on the real slides with
`scripts/senior/build_tissuumaps_project.py` and on the shipped pyramids under
`/data/tissuumaps/channels/`. Donor 6436 (59-plex, 29760 x 47520, 1,414,195,200 px per channel) is
the worked example; the same estimator runs on every donor.

**Marker intensity is not a biological readout in this cohort.** Nothing here is a positivity
threshold, nothing downstream reads any of these numbers, and none of it licenses a statement about
a marker. It is a display stretch, and the question asked of it is only whether it renders what is
in the file.

---

## 1. Where the stretch lives, settled

The shipped `.tmap` files carry `Brightness 0, Gamma 1` on all 59 layers, differing only in `Color`.
That is **by design, not an omission** — `FILTER_DEFAULTS` in the builder says so: the exposure
decision is made during extraction on the full-resolution 16-bit page, where there is far more
signal than in the 8-bit proxies, and is baked into the `.tif` pyramid.

It was applied, and it is the current version. Measured two independent ways:

* Each donor's `channels/<case>/extraction.json` records, for all 59 channels,
  `{"stretch": "robust", "version": 4, "background_k": 3.0, "signal_floor_k": 20.0,
  "signal_quantile": 99.0, "compression": "deflate"}`.
* Stronger, because a manifest is only a claim: for **all 59 channels of 6436** the 8-bit histogram
  predicted by running the v4 robust LUT over the 16-bit source histogram is **bit-identical** to the
  shipped tile's own histogram (largest 256-bin L1 difference 0.0000%). The shipped pixels *are*
  what that code produced — and the same check on the rebuilt donor passes 59 of 59 against the v5
  LUT (section 5), so the method that establishes this is itself confirmed on written pixels.

So "the stretch was never applied" and "the tiles predate it" are both false. The defect is one
level up: **the estimator that picks the cut points never ran as documented on a real slide.**

## 2. The defect

`channel_stats(..., stretch="robust")` is documented to anchor the black point to
`background mode + 3 sigma`, with sigma from a one-sided MAD *below* the mode. On these slides it
does neither.

PhenoCycler subtracts a per-cycle blank and clamps the result at zero, so bin 0 collects every
pixel whose corrected value came out at or below zero. It is a censored pile, not a sample of the
background at level 0 — the mirror image of the saturation bin at the top, which the code already
excludes for exactly that reason. It is also enormous, and it is adjacent to nothing:

| | 6436, all 59 channels |
|---|---|
| share of the frame in bin 0 | **33.2% (Bcl-2) to 91.7% (CD57)** |
| ratio bin 0 : bin 1 | 109x (INS), 444x (CD3e), 1100x (DAPI) |

`np.argmax(counts[:sat_bin])` therefore returned **mode = 0 on 59 channels of 59**. The MAD window
`counts[:mode+1]` was then a single bin, returned 0, and sigma fell to its floor — **sigma = 1.0 on
59 of 59**. At that point nothing is a measurement:

* black point `= mode + 3 sigma` = **3 counts on every channel of every donor**,
* white-point floor `= mode + 20 sigma` = 20 counts,
* blow-out reference `= mode + 40 sigma` = 40 counts,

for a channel whose real background peak sits at 1 (INS) and for one whose real background peak
sits at 119 (CD209) alike.

### Why the test suite never saw it

`tissuumaps_smoke_test.py::make_slide` builds its background as
`rng.normal(200, 60) + rng.gamma(1.2, 70)`. Measured on that fixture: bin 0 holds **0.003–0.006%**
of the frame and the background mode lands at 256–269 — a clean interior peak. The fixture has no
censored pile, so the estimator finds a real mode on it and passes. The one thing that breaks it on
every real slide is the one thing the fixture does not contain.

### What it looked like in the viewer

Two failures at once, both measured as output grey levels 0–255 on 6436:

* **Background lifted into the light.** The median pixel of the whole section rendered at 39
  (CD209), 33 (CD3e), 29 (TP63), 25 (PD-L1), 18 (GranzymeB). At the 90th percentile: Bcl-2 168,
  GranzymeB 164, CD4 147, CD3e 143, TP63 143. A tenth of a pancreas section at or above half
  brightness on a T-cell marker.
* **Signal crushed into flat white.** The 99.9th-percentile pixel rendered at a flat 255 on 12 of 13
  channels checked, with 0.6% of the entire frame saturated on CD3e, ICOS, Bcl-2, TP63, CD209,
  PD-L1 and DAPI.

The builder's own blow-out guard could not see the second one, because it measures against
`mode + 40 sigma` and sigma was the floor. Re-evaluating the same shipped renderings with a working
background model: **20 of 59 channels clip more than 5% of their own signal, and nine clip 100% of
it** — FOXP3, GranzymeB, PD-1, ICOS, PD-L1, CD3e, CD56, CD209, TP63. The guard as it ran flagged
**1 of 59** and reported a median of 1.2%.

## 3. The fix

One change, in `channel_stats`: search the background mode over the bins that sample the
background, which excludes the clamped bottom bin as well as the saturated top bin, and measure the
MAD over the same window.

```python
base_bin = 1 if counts.size > 2 else 0
mode = base_bin + int(np.argmax(counts[base_bin:sat_bin])) if sat_bin > base_bin else base_bin
below = counts[base_bin: mode + 1]
```

Nothing is discarded: bin-0 pixels still map through the LUT to output 0, which is where they were
going anyway. `EXTRACTION_VERSION` goes 4 -> 5, so every cached channel is recognised as stale.

**A second estimator was tried and rejected on measurement.** Where the peak lands on `base_bin`
there is no lower side to measure and sigma keeps its floor. The obvious substitute is the
half-width at half-maximum on the *upper* side. Across the 22 channels of 6436 where both can be
computed it runs a median **1.62x** and up to **14.3x** the lower-side MAD, because the background's
upper tail is the autofluorescence the signal floor exists to step over. It measures the tail, not
the peak, so it is not a substitute — and adopting it would have been choosing a darker picture,
not a better measurement. It is not in the code.

Consequence of that decision: on 6436, **34 of 59** channels keep sigma = 1.0 and render as before.
All 34 already rendered their median pixel at output 0, and only SOX2 put its 90th percentile above
output 100. The **25** channels that gain a measured sigma are exactly the ones whose background had
been lifted into the light.

## 4. Before and after, per channel (6436)

`black`/`white` are 16-bit cut points; `median px`/`p90 px` are the 0–255 output level the section's
median and 90th-percentile pixel render at; `% frame >=128` is how much of the frame renders at or
above half brightness; `% signal white` is the fraction of unambiguously-bright pixels rendered flat
white, both columns evaluated against the corrected background model so they are comparable.

| channel | bin-0 % | bg mode | sigma | black | white | median px | p90 px | % frame >=128 | % signal white |
|---|---|---|---|---|---|---|---|---|---|
| DAPI | 33.2 | 0 -> 27 | 1 -> 10 | 3 -> 58 | 1601 -> 1841 | 6 -> 0 | 115 -> 95 | 8.51 -> 6.08 | 4.0 -> 1.5 |
| IAPP | 87.5 | 0 -> 1 | 1 -> 1 | 3 -> 4 | 215 -> 216 | 0 -> 0 | 7 -> 6 | 0.94 -> 0.92 | 1.6 -> 1.6 |
| CD107a | 51.4 | 0 -> 1 | 1 -> 1 | 3 -> 4 | 752 -> 753 | 0 -> 0 | 79 -> 79 | 3.89 -> 3.87 | 1.1 -> 1.1 |
| CD4 | 40.2 | 0 -> 13 | 1 -> 9 | 3 -> 39 | 763 -> 823 | 7 -> 0 | 147 -> 130 | 16.47 -> 10.73 | 2.8 -> 1.9 |
| CD31 | 62.2 | 0 -> 1 | 1 -> 1 | 3 -> 4 | 4550 -> 4598 | 0 -> 0 | 5 -> 5 | 0.62 -> 0.61 | 1.4 -> 1.3 |
| SMA | 67.3 | 0 -> 1 | 1 -> 1 | 3 -> 4 | 6503 -> 6673 | 0 -> 0 | 1 -> 1 | 0.23 -> 0.22 | 2.5 -> 2.3 |
| CD68 | 48.1 | 0 -> 4 | 1 -> 2 | 3 -> 8 | 611 -> 657 | 0 -> 0 | 61 -> 55 | 1.36 -> 1.16 | 1.4 -> 1.2 |
| CD44 | 45.3 | 0 -> 1 | 1 -> 1 | 3 -> 4 | 12833 -> 12851 | 0 -> 0 | 54 -> 54 | 2.95 -> 2.94 | 1.1 -> 1.1 |
| Vimentin | 45.3 | 0 -> 1 | 1 -> 1 | 3 -> 4 | 6604 -> 6618 | 0 -> 0 | 27 -> 27 | 1.81 -> 1.80 | 1.1 -> 1.1 |
| CD99 | 45.9 | 0 -> 2 | 1 -> 1 | 3 -> 5 | 2452 -> 2485 | 0 -> 0 | 23 -> 22 | 0.81 -> 0.81 | 1.2 -> 1.2 |
| ECAD | 56.0 | 0 -> 1 | 1 -> 1 | 3 -> 4 | 5311 -> 5321 | 0 -> 0 | 42 -> 42 | 2.20 -> 2.20 | 1.1 -> 1.1 |
| IDO1 | 73.4 | 0 -> 1 | 1 -> 1 | 3 -> 4 | 108 -> 109 | 0 -> 0 | 44 -> 41 | 1.19 -> 1.12 | 3.5 -> 3.4 |
| CD66 | 44.0 | 0 -> 16 | 1 -> 10 | 3 -> 47 | 4759 -> 13222 | 1 -> 0 | 9 -> 2 | 0.92 -> 0.27 | 21.1 -> 2.2 |
| INS | 85.2 | 0 -> 1 | 1 -> 1 | 3 -> 4 | 330 -> 400 | 0 -> 0 | 3 -> 2 | 0.05 -> 0.05 | 5.1 -> 4.7 |
| Ker8-18 | 46.1 | 0 -> 1 | 1 -> 1 | 3 -> 4 | 4613 -> 4620 | 0 -> 0 | 64 -> 64 | 3.51 -> 3.50 | 1.1 -> 1.1 |
| CD11c | 51.4 | 0 -> 8 | 1 -> 4 | 3 -> 21 | 470 -> 1347 | 0 -> 0 | 35 -> 9 | 0.90 -> 0.18 | 23.7 -> 3.6 |
| CD38 | 80.0 | 0 -> 1 | 1 -> 1 | 3 -> 4 | 339 -> 340 | 0 -> 0 | 33 -> 32 | 2.15 -> 2.11 | 1.3 -> 1.3 |
| HLADR | 70.6 | 0 -> 1 | 1 -> 1 | 3 -> 4 | 164 -> 169 | 0 -> 0 | 40 -> 37 | 0.51 -> 0.47 | 3.1 -> 2.9 |
| CD34 | 36.7 | 0 -> 23 | 1 -> 15 | 3 -> 67 | 6193 -> 9920 | 2 -> 0 | 23 -> 13 | 1.95 -> 0.88 | 6.1 -> 1.9 |
| iNOS | 79.5 | 0 -> 1 | 1 -> 1 | 3 -> 4 | 202 -> 203 | 0 -> 0 | 35 -> 33 | 1.53 -> 1.48 | 1.6 -> 1.6 |
| M2Gal3 | 55.1 | 0 -> 1 | 1 -> 1 | 3 -> 4 | 379 -> 381 | 0 -> 0 | 55 -> 54 | 1.60 -> 1.57 | 1.3 -> 1.3 |
| B3TUBB | 47.5 | 0 -> 11 | 1 -> 7 | 3 -> 33 | 351 -> 4525 | 1 -> 0 | 74 -> 4 | 2.19 -> 0.08 | 85.7 -> 7.0 |
| CD8a | 75.6 | 0 -> 1 | 1 -> 1 | 3 -> 4 | 237 -> 249 | 0 -> 0 | 41 -> 39 | 0.76 -> 0.61 | 1.6 -> 1.6 |
| PCNA | 55.6 | 0 -> 1 | 1 -> 1 | 3 -> 4 | 507 -> 508 | 0 -> 0 | 68 -> 68 | 3.48 -> 3.45 | 1.1 -> 1.1 |
| FOXP3 | 41.1 | 0 -> 23 | 1 -> 15 | 3 -> 67 | 441 -> 591 | 12 -> 0 | 137 -> 84 | 11.98 -> 3.01 | 100.0 -> 100.0 |
| CD57 | 91.7 | 0 -> 1 | 1 -> 1 | 3 -> 4 | 88 -> 94 | 0 -> 0 | 0 -> 0 | 0.22 -> 0.15 | 6.5 -> 6.0 |
| Ki67 | 67.9 | 0 -> 1 | 1 -> 1 | 3 -> 4 | 107 -> 109 | 0 -> 0 | 71 -> 68 | 2.67 -> 2.50 | 2.6 -> 2.5 |
| GranzymeB | 34.4 | 0 -> 43 | 1 -> 19 | 3 -> 100 | 735 -> 889 | 18 -> 0 | 164 -> 121 | 21.24 -> 8.21 | 100.0 -> 51.8 |
| HLA-A | 65.7 | 0 -> 1 | 1 -> 1 | 3 -> 4 | 1439 -> 1443 | 0 -> 0 | 37 -> 37 | 1.95 -> 1.93 | 1.2 -> 1.2 |
| CD20 | 57.8 | 0 -> 1 | 1 -> 1 | 3 -> 4 | 430 -> 432 | 0 -> 0 | 52 -> 51 | 1.54 -> 1.52 | 1.3 -> 1.3 |
| CollagenIV | 34.4 | 0 -> 47 | 1 -> 21 | 3 -> 109 | 2014 -> 3286 | 7 -> 0 | 75 -> 39 | 3.64 -> 1.00 | 12.8 -> 3.0 |
| VISTA | 64.0 | 0 -> 1 | 1 -> 1 | 3 -> 4 | 1596 -> 1604 | 0 -> 0 | 23 -> 22 | 1.04 -> 1.03 | 1.2 -> 1.2 |
| Pan_Cytokeratin | 55.9 | 0 -> 1 | 1 -> 1 | 3 -> 4 | 3147 -> 3151 | 0 -> 0 | 43 -> 43 | 2.51 -> 2.51 | 1.1 -> 1.1 |
| PD-1 | 49.3 | 0 -> 14 | 1 -> 9 | 3 -> 40 | 299 -> 524 | 0 -> 0 | 87 -> 34 | 4.12 -> 0.49 | 100.0 -> 14.3 |
| LAG3 | 64.2 | 0 -> 1 | 1 -> 1 | 3 -> 4 | 338 -> 339 | 0 -> 0 | 84 -> 84 | 4.93 -> 4.86 | 1.3 -> 1.2 |
| SST | 69.5 | 0 -> 1 | 1 -> 1 | 3 -> 4 | 300 -> 312 | 0 -> 0 | 31 -> 29 | 0.41 -> 0.38 | 2.0 -> 1.9 |
| TCF-1 | 39.2 | 0 -> 20 | 1 -> 12 | 3 -> 55 | 1068 -> 1339 | 5 -> 0 | 63 -> 42 | 3.34 -> 2.10 | 13.6 -> 2.9 |
| TOX | 52.8 | 0 -> 3 | 1 -> 2 | 3 -> 7 | 464 -> 475 | 0 -> 0 | 108 -> 105 | 6.78 -> 6.24 | 1.4 -> 1.2 |
| Caveolin | 62.2 | 0 -> 1 | 1 -> 1 | 3 -> 4 | 3205 -> 3226 | 0 -> 0 | 6 -> 6 | 1.09 -> 1.08 | 1.5 -> 1.5 |
| ICOS | 35.2 | 0 -> 44 | 1 -> 19 | 3 -> 101 | 789 -> 1424 | 14 -> 0 | 102 -> 41 | 6.07 -> 0.68 | 100.0 -> 8.8 |
| CD163 | 82.2 | 0 -> 1 | 1 -> 1 | 3 -> 4 | 195 -> 198 | 0 -> 0 | 21 -> 20 | 0.64 -> 0.60 | 2.0 -> 1.9 |
| EpCAM | 53.7 | 0 -> 1 | 1 -> 1 | 3 -> 4 | 898 -> 898 | 0 -> 0 | 88 -> 88 | 4.68 -> 4.68 | 1.0 -> 1.0 |
| PD-L1 | 37.0 | 0 -> 58 | 1 -> 27 | 3 -> 138 | 453 -> 1801 | 25 -> 0 | 102 -> 7 | 6.02 -> 0.03 | 100.0 -> 14.8 |
| CD79a | 54.7 | 0 -> 5 | 1 -> 3 | 3 -> 13 | 319 -> 339 | 0 -> 0 | 118 -> 106 | 8.41 -> 6.65 | 2.7 -> 1.8 |
| Keratin 5 | 80.1 | 0 -> 1 | 1 -> 1 | 3 -> 4 | 121 -> 122 | 0 -> 0 | 32 -> 30 | 1.35 -> 1.28 | 2.3 -> 2.3 |
| CD3e | 35.5 | 0 -> 47 | 1 -> 22 | 3 -> 113 | 339 -> 4485 | 33 -> 0 | 143 -> 5 | 13.75 -> 0.02 | 100.0 -> 1.7 |
| GCG | 52.8 | 0 -> 2 | 1 -> 1 | 3 -> 5 | 3950 -> 4024 | 0 -> 0 | 6 -> 6 | 0.79 -> 0.78 | 1.4 -> 1.4 |
| Beta-actin | 47.9 | 0 -> 1 | 1 -> 1 | 3 -> 4 | 2622 -> 2625 | 0 -> 0 | 59 -> 59 | 2.71 -> 2.70 | 1.1 -> 1.1 |
| Bcl-2 | 33.2 | 0 -> 109 | 1 -> 28 | 3 -> 193 | 2345 -> 2527 | 15 -> 0 | 168 -> 148 | 24.01 -> 16.48 | 3.1 -> 1.7 |
| MPO | 50.5 | 0 -> 4 | 1 -> 2 | 3 -> 8 | 2543 -> 2755 | 0 -> 0 | 14 -> 13 | 0.92 -> 0.85 | 1.5 -> 1.3 |
| CD56 | 37.3 | 0 -> 42 | 1 -> 22 | 3 -> 108 | 670 -> 2652 | 12 -> 0 | 100 -> 16 | 6.77 -> 0.09 | 100.0 -> 11.8 |
| CD39 | 81.6 | 0 -> 1 | 1 -> 1 | 3 -> 4 | 726 -> 738 | 0 -> 0 | 6 -> 5 | 0.34 -> 0.33 | 2.0 -> 2.0 |
| Iba1 | 35.6 | 0 -> 33 | 1 -> 18 | 3 -> 86 | 1545 -> 4989 | 5 -> 0 | 26 -> 4 | 1.42 -> 0.22 | 39.0 -> 1.7 |
| SOX2 | 54.2 | 0 -> 1 | 1 -> 1 | 3 -> 4 | 515 -> 516 | 0 -> 0 | 106 -> 106 | 6.41 -> 6.34 | 1.2 -> 1.2 |
| CD209 | 33.9 | 0 -> 119 | 1 -> 31 | 3 -> 212 | 749 -> 4313 | 39 -> 0 | 131 -> 11 | 10.44 -> 0.03 | 100.0 -> 11.2 |
| PDPN | 65.0 | 0 -> 1 | 1 -> 1 | 3 -> 4 | 381 -> 384 | 0 -> 0 | 44 -> 43 | 0.93 -> 0.91 | 1.5 -> 1.4 |
| CD11b | 34.6 | 0 -> 28 | 1 -> 15 | 3 -> 72 | 1550 -> 3087 | 5 -> 0 | 39 -> 14 | 2.41 -> 0.54 | 16.8 -> 2.1 |
| CD206 | 39.6 | 0 -> 20 | 1 -> 13 | 3 -> 60 | 909 -> 2377 | 8 -> 0 | 103 -> 34 | 3.95 -> 0.40 | 32.6 -> 12.3 |
| TP63 | 34.2 | 0 -> 42 | 1 -> 19 | 3 -> 99 | 379 -> 2392 | 29 -> 0 | 143 -> 13 | 14.83 -> 0.02 | 100.0 -> 6.3 |

Summary over the 59 channels:

| | shipped | rebuilt |
|---|---|---|
| background mode | 0 on 59/59 | 1–119, measured |
| sigma | 1.0 on 59/59 | 1.0–31.1 (25 measured, 34 at the floor) |
| black point | 3 on 59/59 | 4–212 |
| median pixel output | up to 39 | **0 on all 59** |
| mean % of frame at >=128 | 4.15% | **2.17%** |
| channels clipping >5% of their own signal | 20 | **11** |
| channels clipping 100% of their own signal | 9 | 1 (FOXP3) |

Verified on rebuilt pixels, 6436 CD3e (the rebuilt tile is bit-identical to the predicted LUT):

| | shipped tile | rebuilt tile |
|---|---|---|
| pixels at output 0 | 35.78% | 74.61% |
| pixels at output 255 | 0.6178% | 0.0017% |
| mean output | 52.65 | 1.23 |
| median / p90 output | 33 / 143 | 0 / 5 |
| % of frame >=128 | 13.95% | 0.020% |

---

## 5. Before and after on REBUILT pixels (6539)

Section 4 compares two renderings of 6436 computed from the same 16-bit source histogram, with the
shipped side proven bit-identical to the tile on disk. This section is the same comparison on the
first donor actually rebuilt, so the "after" column is **measured on the pixels the viewer now
serves**, not predicted.

Donor 6539, 59 channels, 24960 x 23760 = 593,049,600 px per channel. One histogram pass over the
16-bit source page and one over the rebuilt 8-bit pyramid level 0, per channel.

**Verification first.** For each channel the 8-bit histogram predicted by running today's LUT over
the source histogram was compared against the rebuilt tile's own histogram:
**59 of 59 are bit-identical** (largest 256-bin L1 difference 0.000000%). The pixels on disk are
what this code produces, so every column below is a measurement of the shipped rendering.

Columns as in section 4: `black`/`white` are 16-bit cut points; `median px`/`p90 px` are the 0-255
output level the section's median and 90th-percentile pixel render at; `% signal white` is the
fraction of unambiguously-bright pixels (source level above mode + 40 sigma, from the working
background model, so both columns score the same population) rendered flat white.

| channel | bin-0 % | bg mode | sigma | black | white | median px | p90 px | % frame >=128 | % signal white |
|---|---|---|---|---|---|---|---|---|---|
| DAPI | 30.6 | 0 -> 27 | 1 -> 10 | 3 -> 58 | 1604 -> 1994 | 5 -> 0 | 59 -> 41 | 5.04 -> 3.32 | 7.1 -> 1.7 |
| IAPP | 90.3 | 0 -> 1 | 1 -> 1 | 3 -> 4 | 222 -> 224 | 0 -> 0 | 0 -> 0 | 0.55 -> 0.54 | 1.8 -> 1.8 |
| CD107a | 56.2 | 0 -> 1 | 1 -> 1 | 3 -> 4 | 847 -> 849 | 0 -> 0 | 58 -> 58 | 2.01 -> 1.99 | 1.1 -> 1.1 |
| CD4 | 36.0 | 0 -> 21 | 1 -> 13 | 3 -> 61 | 1139 -> 1393 | 6 -> 0 | 151 -> 118 | 17.23 -> 7.18 | 3.0 -> 1.6 |
| CD31 | 56.7 | 0 -> 1 | 1 -> 1 | 3 -> 4 | 3404 -> 3445 | 0 -> 0 | 10 -> 10 | 0.65 -> 0.65 | 1.3 -> 1.3 |
| SMA | 53.5 | 0 -> 3 | 1 -> 2 | 3 -> 7 | 8872 -> 9746 | 0 -> 0 | 1 -> 1 | 1.33 -> 1.19 | 3.1 -> 1.9 |
| CD68 | 46.5 | 0 -> 11 | 1 -> 7 | 3 -> 33 | 574 -> 3483 | 1 -> 0 | 38 -> 4 | 0.85 -> 0.08 | 47.9 -> 3.3 |
| CD44 | 55.8 | 0 -> 1 | 1 -> 1 | 3 -> 4 | 5501 -> 5516 | 0 -> 0 | 34 -> 34 | 1.86 -> 1.85 | 1.1 -> 1.1 |
| Vimentin | 43.9 | 0 -> 1 | 1 -> 1 | 3 -> 4 | 10056 -> 10083 | 0 -> 0 | 17 -> 17 | 1.66 -> 1.66 | 1.1 -> 1.1 |
| CD99 | 37.0 | 0 -> 12 | 1 -> 7 | 3 -> 34 | 1601 -> 3404 | 2 -> 0 | 40 -> 17 | 0.95 -> 0.41 | 5.9 -> 2.4 |
| ECAD | 61.3 | 0 -> 1 | 1 -> 1 | 3 -> 4 | 6415 -> 6429 | 0 -> 0 | 38 -> 38 | 1.99 -> 1.98 | 1.1 -> 1.1 |
| IDO1 | 64.0 | 0 -> 1 | 1 -> 1 | 3 -> 4 | 181 -> 182 | 0 -> 0 | 54 -> 53 | 1.43 -> 1.39 | 2.1 -> 2.1 |
| CD66 | 38.2 | 0 -> 26 | 1 -> 16 | 3 -> 74 | 306 -> 5218 | 15 -> 0 | 129 -> 4 | 10.24 -> 0.01 | 100.0 -> 5.0 |
| INS | 67.4 | 0 -> 1 | 1 -> 1 | 3 -> 4 | 6701 -> 6914 | 0 -> 0 | 1 -> 1 | 0.29 -> 0.28 | 2.8 -> 2.7 |
| Ker8-18 | 43.6 | 0 -> 4 | 1 -> 2 | 3 -> 8 | 7075 -> 7398 | 0 -> 0 | 41 -> 39 | 1.99 -> 1.84 | 1.4 -> 1.2 |
| CD11c | 34.0 | 0 -> 35 | 1 -> 18 | 3 -> 88 | 451 -> 2226 | 15 -> 0 | 63 -> 3 | 1.98 -> 0.07 | 100.0 -> 3.7 |
| CD38 | 77.5 | 0 -> 1 | 1 -> 1 | 3 -> 4 | 376 -> 378 | 0 -> 0 | 22 -> 21 | 1.41 -> 1.39 | 1.6 -> 1.5 |
| HLADR | 63.3 | 0 -> 1 | 1 -> 1 | 3 -> 4 | 281 -> 285 | 0 -> 0 | 43 -> 42 | 0.74 -> 0.71 | 1.8 -> 1.7 |
| CD34 | 32.1 | 0 -> 34 | 1 -> 18 | 3 -> 87 | 7155 -> 12196 | 2 -> 0 | 23 -> 12 | 2.22 -> 0.87 | 7.0 -> 1.8 |
| iNOS | 71.0 | 0 -> 1 | 1 -> 1 | 3 -> 4 | 229 -> 231 | 0 -> 0 | 58 -> 56 | 2.69 -> 2.56 | 1.5 -> 1.5 |
| M2Gal3 | 47.9 | 0 -> 1 | 1 -> 1 | 3 -> 4 | 990 -> 993 | 0 -> 0 | 69 -> 69 | 1.82 -> 1.81 | 1.1 -> 1.1 |
| B3TUBB | 37.0 | 0 -> 24 | 1 -> 15 | 3 -> 68 | 596 -> 16407 | 8 -> 0 | 56 -> 1 | 1.03 -> 0.06 | 100.0 -> 1.8 |
| CD8a | 68.4 | 0 -> 1 | 1 -> 1 | 3 -> 4 | 309 -> 319 | 0 -> 0 | 65 -> 62 | 1.36 -> 1.14 | 1.3 -> 1.3 |
| PCNA | 55.5 | 0 -> 1 | 1 -> 1 | 3 -> 4 | 981 -> 983 | 0 -> 0 | 52 -> 52 | 1.94 -> 1.92 | 1.1 -> 1.1 |
| FOXP3 | 34.2 | 0 -> 37 | 1 -> 21 | 3 -> 99 | 526 -> 746 | 16 -> 0 | 146 -> 80 | 14.49 -> 2.24 | 100.0 -> 100.0 |
| CD57 | 79.0 | 0 -> 1 | 1 -> 1 | 3 -> 4 | 77 -> 78 | 0 -> 0 | 52 -> 48 | 2.57 -> 2.39 | 3.9 -> 3.6 |
| Ki67 | 57.9 | 0 -> 1 | 1 -> 1 | 3 -> 4 | 119 -> 119 | 0 -> 0 | 114 -> 113 | 8.20 -> 7.86 | 1.6 -> 1.6 |
| GranzymeB | 30.8 | 0 -> 74 | 1 -> 24 | 3 -> 145 | 1115 -> 1276 | 19 -> 0 | 169 -> 134 | 19.72 -> 11.54 | 49.0 -> 14.3 |
| HLA-A | 61.0 | 0 -> 1 | 1 -> 1 | 3 -> 4 | 1416 -> 1425 | 0 -> 0 | 20 -> 19 | 1.01 -> 1.00 | 1.4 -> 1.3 |
| CD20 | 52.3 | 0 -> 1 | 1 -> 1 | 3 -> 4 | 442 -> 444 | 0 -> 0 | 76 -> 75 | 2.48 -> 2.45 | 1.2 -> 1.1 |
| CollagenIV | 31.2 | 0 -> 63 | 1 -> 27 | 3 -> 143 | 1665 -> 2886 | 10 -> 0 | 90 -> 42 | 4.73 -> 0.91 | 30.7 -> 4.5 |
| VISTA | 54.2 | 0 -> 1 | 1 -> 1 | 3 -> 4 | 2191 -> 2203 | 0 -> 0 | 29 -> 29 | 1.16 -> 1.15 | 1.2 -> 1.2 |
| Pan_Cytokeratin | 55.7 | 0 -> 1 | 1 -> 1 | 3 -> 4 | 5247 -> 5260 | 0 -> 0 | 31 -> 31 | 1.75 -> 1.74 | 1.1 -> 1.1 |
| PD-1 | 34.1 | 0 -> 45 | 1 -> 22 | 3 -> 111 | 335 -> 1850 | 28 -> 0 | 119 -> 7 | 8.49 -> 0.01 | 100.0 -> 15.9 |
| LAG3 | 61.9 | 0 -> 1 | 1 -> 1 | 3 -> 4 | 493 -> 494 | 0 -> 0 | 72 -> 71 | 3.58 -> 3.55 | 1.2 -> 1.2 |
| SST | 64.1 | 0 -> 1 | 1 -> 1 | 3 -> 4 | 215 -> 218 | 0 -> 0 | 52 -> 50 | 1.33 -> 1.25 | 1.8 -> 1.8 |
| TCF-1 | 32.9 | 0 -> 37 | 1 -> 19 | 3 -> 94 | 1741 -> 2686 | 5 -> 0 | 52 -> 26 | 2.30 -> 1.08 | 23.8 -> 2.6 |
| TOX | 47.4 | 0 -> 9 | 1 -> 6 | 3 -> 26 | 532 -> 585 | 1 -> 0 | 124 -> 107 | 9.40 -> 6.26 | 3.7 -> 2.1 |
| Caveolin | 53.3 | 0 -> 1 | 1 -> 1 | 3 -> 4 | 5206 -> 5226 | 0 -> 0 | 7 -> 7 | 1.30 -> 1.29 | 1.2 -> 1.2 |
| ICOS | 31.4 | 0 -> 61 | 1 -> 25 | 3 -> 136 | 1023 -> 2342 | 14 -> 0 | 77 -> 20 | 3.67 -> 0.40 | 100.0 -> 4.7 |
| CD163 | 69.0 | 0 -> 1 | 1 -> 1 | 3 -> 4 | 174 -> 176 | 0 -> 0 | 57 -> 55 | 2.19 -> 2.12 | 1.9 -> 1.8 |
| EpCAM | 44.5 | 0 -> 4 | 1 -> 2 | 3 -> 8 | 1770 -> 1815 | 1 -> 0 | 87 -> 85 | 4.88 -> 4.57 | 1.3 -> 1.1 |
| PD-L1 | 31.7 | 0 -> 96 | 1 -> 33 | 3 -> 193 | 460 -> 2028 | 45 -> 0 | 107 -> 0 | 6.58 -> 0.01 | 100.0 -> 22.8 |
| CD79a | 42.0 | 0 -> 22 | 1 -> 15 | 3 -> 66 | 433 -> 1667 | 9 -> 0 | 134 -> 26 | 11.26 -> 0.08 | 100.0 -> 27.7 |
| Keratin 5 | 61.5 | 0 -> 1 | 1 -> 1 | 3 -> 4 | 163 -> 164 | 0 -> 0 | 59 -> 57 | 2.57 -> 2.48 | 2.0 -> 2.0 |
| CD3e | 31.3 | 0 -> 63 | 1 -> 24 | 3 -> 134 | 424 -> 10191 | 35 -> 0 | 143 -> 3 | 14.47 -> 0.03 | 100.0 -> 1.6 |
| GCG | 51.9 | 0 -> 4 | 1 -> 2 | 3 -> 8 | 417 -> 578 | 0 -> 0 | 46 -> 31 | 0.93 -> 0.45 | 2.2 -> 1.7 |
| Beta-actin | 41.1 | 0 -> 6 | 1 -> 3 | 3 -> 14 | 2742 -> 2913 | 1 -> 0 | 56 -> 52 | 3.02 -> 2.64 | 1.4 -> 1.2 |
| Bcl-2 | 30.6 | 0 -> 139 | 1 -> 38 | 3 -> 254 | 3392 -> 3811 | 12 -> 0 | 168 -> 142 | 21.57 -> 14.52 | 3.2 -> 1.4 |
| MPO | 47.0 | 0 -> 11 | 1 -> 7 | 3 -> 33 | 268 -> 529 | 2 -> 0 | 132 -> 55 | 10.83 -> 0.29 | 100.0 -> 42.9 |
| CD56 | 31.5 | 0 -> 68 | 1 -> 24 | 3 -> 139 | 704 -> 4419 | 22 -> 0 | 71 -> 3 | 5.03 -> 0.10 | 100.0 -> 5.8 |
| CD39 | 48.9 | 0 -> 29 | 1 -> 19 | 3 -> 86 | 929 -> 2763 | 0 -> 0 | 20 -> 0 | 1.31 -> 0.14 | 74.0 -> 2.9 |
| Iba1 | 30.8 | 0 -> 71 | 1 -> 24 | 3 -> 142 | 1331 -> 6627 | 13 -> 0 | 38 -> 2 | 1.17 -> 0.14 | 78.9 -> 1.5 |
| SOX2 | 45.5 | 0 -> 16 | 1 -> 10 | 3 -> 47 | 640 -> 779 | 3 -> 0 | 114 -> 84 | 7.77 -> 3.63 | 14.5 -> 4.9 |
| CD209 | 30.7 | 0 -> 166 | 1 -> 40 | 3 -> 286 | 1079 -> 17138 | 37 -> 0 | 96 -> 2 | 5.49 -> 0.05 | 100.0 -> 2.4 |
| PDPN | 45.4 | 0 -> 21 | 1 -> 15 | 3 -> 65 | 306 -> 3784 | 7 -> 0 | 117 -> 5 | 7.89 -> 0.02 | 100.0 -> 2.4 |
| CD11b | 30.6 | 0 -> 57 | 1 -> 18 | 3 -> 110 | 5935 -> 8256 | 2 -> 0 | 45 -> 29 | 3.61 -> 1.86 | 5.5 -> 1.6 |
| CD206 | 35.8 | 0 -> 24 | 1 -> 15 | 3 -> 68 | 2065 -> 6142 | 4 -> 0 | 60 -> 18 | 0.91 -> 0.43 | 15.3 -> 6.6 |

Summary over the 59 channels:

| TP63 | 30.7 | 0 -> 60 | 1 -> 21 | 3 -> 122 | 442 -> 1390 | 34 -> 0 | 177 -> 37 | 22.61 -> 0.03 | 100.0 -> 24.9 |

| | shipped (v4) | rebuilt (v5) |
|---|---|---|
| background mode | 0–0 (59/59 at 0) | 1–166 (0/59 at 0) |
| sigma at its 1.0 floor | 59/59 | 25/59 |
| black point | 3–3 | 4–286 |
| median pixel output > 0 | 30/59 | 0/59 |
| mean % of frame at >=128 | 4.81% | 1.93% |
| channels clipping >5% of their own signal | 26 | 10 |
| channels clipping 100% of their own signal | 14 | 1 |

The shape of the defect is the same as on 6436 and somewhat worse: **30 of 59** channels rendered
the median pixel of the whole section above output 0 (up to 45, on PD-L1), and **14 of 59** rendered
100% of their own unambiguous signal as flat white. After the rebuild the median pixel renders at 0
on all 59, and one channel — FOXP3 — still clips 100%.

## 6. The docstring's comparison, reproduced on real channels

`channel_stats()` documents the robust stretch against a fixed `percentile 0.5/99.5` on four
synthetic cases: dense, sparse, high background, saturated. Both estimators were run over the same
6539 source histograms, and each rendering scored against the same "certainly signal" population.

### percentile 0.5/99.5 vs robust, on real channels

(no saturated counterpart: the top histogram bin holds 0.0000% of the frame at most, over all 59 channels)

| channel (role) | bin-0 % | detector-sat % | percentile black/white | robust black/white | median px pct -> rob | % frame >=128 pct -> rob | % signal white pct -> rob |
|---|---|---|---|---|---|---|---|
| M2Gal3 (dense) | 47.9 | 0.000 | 0/911 | 4/993 | 1 -> 0 | 2.25 -> 1.81 | 1.4 -> 1.1 |
| CD39 (sparse) | 48.9 | 0.000 | 0/827 | 86/2763 | 1 -> 0 | 1.54 -> 0.14 | 93.6 -> 2.9 |
| CD209 (high background) | 30.7 | 0.000 | 0/1336 | 286/17138 | 31 -> 0 | 3.16 -> 0.05 | 100.0 -> 2.4 |

Over all channels:
* percentile: median pixel renders above 0 on 32/59, mean % of frame >=128 5.23%, 32 channels clip >5% of their own signal (13 clip 100%).

The fourth synthetic case does not occur here. `hist_find` sizes the histogram to the data range, so
its top bin is the channel's own maximum rather than a 16-bit ceiling, and it holds **0.0000% of the
frame on all 59 channels of 6539 and all 59 of 6436** — nothing on these slides piles up against
detector saturation. The saturation-bin exclusion in `channel_stats()` is therefore doing no work on
this cohort; it is the *bottom* bin, the clamped one, that carries 30-90% of every frame. That
asymmetry is the whole defect in section 2, and it is the reverse of what the code was written
expecting.

`percentile` is not a fallback to reach for: on 32 of 59 channels it renders the median pixel of the
section above black, and on 32 of 59 it renders more than 5% of the channel's own signal as flat
white (13 of those at 100%). It is retained only as the comparison `--stretch percentile`.

## 7. What the fix does NOT settle

The blow-out guard now has a valid background model and, for the first time, something to say.
After the rebuild **11 of 59** channels still clip more than 5% of their own signal on 6436:
FOXP3 100%, GranzymeB 51.8%, PD-L1 14.8%, PD-1 14.3%, CD206 12.3%, CD56 11.8%, CD209 11.2%,
ICOS 8.8%, B3TUBB 7.0%, TP63 6.3%, CD57 6.0%. On FOXP3 the white point lands *below* the
"certainly signal" reference, which is the histogram saying that channel's population is not
separable from its own background tail at these settings.

`--signal-floor-k` was **not** raised to clear those. Turning a verdict green by moving the
parameter it is computed from is tuning to the verdict; the guard has been blind since it was
written and its first honest reading belongs in front of the maintainer, not smoothed away.

## 8. Two other things found in the same file

**A live bug in commit 4335647** (`Enhance layer opacity management`, 2026-09-15). In
`build_project()` the nuclear-opacity test read `name.upper()`, where `name` is the loop variable
left over from the layer loop above and therefore holds the **last channel of the panel** (TP63 on a
59-plex), not the marker being switched on. `DEFAULT_ON` is ordered INS, GCG, SST, DAPI, so the one
layer that decision affects — DAPI — was being tested against TP63. Measured: rebuilding 6436's
project with HEAD's code produced a tmap identical to the shipped one **except DAPI at 0.70 instead
of 0.35**, which is the tissue-wide wash those constants were added to remove. Now reads
`channels[idx].upper()`.

**IAPP, and one other.** The marker registry declares two channels `kind: "excluded"` — IAPP with
`uses: ["failed_marker"]`, and CD11b with `uses: ["assay_qc"]`. Both ship as rendered layers on
every donor.

*Decision: ship them, labelled with the registry's own reason, default-off as they already were.*
The layer names now read `IAPP [FAILED STAIN - not readable]` and
`CD11b [ASSAY QC ONLY - not an identity marker]`. Dropping them was considered and rejected: a
59-plex slide that opens with 58 layers is a silent edit, and whoever goes looking for IAPP would
have no way to tell an excluded channel from a build that lost one. The wording is read from
`marker_registry.json` at build time rather than written here, an unknown reason is passed through
verbatim, and a registry that cannot be read leaves every layer unannotated and logs that it did.

This is **display text only**. Channels are still resolved by marker name against each slide's own
OME-XML, never by index, and nothing reads a layer name back — TissUUmaps assigns
`td_name.innerHTML` and keys every element id off the layer index
(`overlayUtils.addLayerSettings`). The hue is still hashed from the bare marker name, so an excluded
channel keeps the same colour it has on every other donor.

Worth a second opinion, not settled here: `core/phenocycler/scripts/typing/lineages.py` names the
two held-back channels as "IAPP and DAPI", and actively uses CD11b — while the registry marks CD11b
`excluded`. The two disagree about CD11b.

## 9. Rebuilding

The channel pyramids have to be rewritten for the change to reach the viewer; the `.tmap` alone
carries none of the stretch.

```bash
cd /home/smith6jt/IO60panc2nd/Islet-Explorer-Senior
.venv/bin/python scripts/senior/build_tissuumaps_project.py \
    --image /data/islet_ome_tiff/<case>.ome.tiff --out-dir /data/tissuumaps
```

`--regions` is deliberately omitted: `copy_regions()` would rewrite
`/data/tissuumaps/regions/<case>.geojson`, which another lane owns and is separately correcting, and
pointing it at the already-scaled copy would scale it twice. The builder then emits an empty
`regionFiles`, so the block is restored verbatim from the pre-rebuild backup immediately afterwards
and the islets keep auto-loading. Pre-rebuild copies of all ten `.tmap` files are in the session
scratchpad.

Cost, measured: about 27 s per channel on an idle machine and about 2 min per channel under load,
so 30 min to 2 h per donor, and roughly 18 GB per donor on disk. The rebuilt files are *smaller*
than the ones they replace (6436 CD3e 446 MB -> 304 MB; 6539 DAPI 346 MB -> 258 MB), because less
background noise is being encoded.

Nothing under `/data/tissuumaps/channels/` is deleted: each channel is written to `<name>.tif.part`
and moved into place.

## 10. Left for the maintainer

* **The smoke-test fixture still cannot catch this class of defect.**
  `tissuumaps_smoke_test.py::make_slide` is outside this lane's ownership and was not edited. Its
  background (`normal(200, 60) + gamma(1.2, 70)`) puts 0.003–0.006% of the frame in bin 0 where the
  real slides put 33–92%, so the estimator finds a clean interior mode on it either way and the
  regression that shipped is invisible to the suite. Adding a clamped, censored floor to the fixture
  — a background whose mean sits at or below zero, so a large share of it piles into bin 0 — would
  make the failure reproducible in the test rather than only on the slides.
* **`lineages.py` and the registry disagree about CD11b.** The registry has it `kind: "excluded"`,
  `uses: ["assay_qc"]`; `lineages.py` names the held-back channels as "IAPP and DAPI" and uses CD11b
  for typing. The viewer now labels CD11b from the registry. Whichever is right, one of them should
  move.
* **Eleven channels still clip more than 5% of their own signal** at the default `--signal-floor-k`
  (section 7). That parameter was deliberately left alone.
