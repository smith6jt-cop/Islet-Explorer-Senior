# UniFORM 22-donor alignment — 2026-09-16

> **CORRECTED 2026-09-18.** "Donors 6457 and 6579 are excluded" was true when this ran and is not
> now: **only 6457 is centrally excluded**. **6579's exclusion was lifted 2026-09-17** and it is
> recorded in `REINSTATED` in `core/phenocycler/phenocycler/cohort.py`; the pancreatitis it was
> excluded for is still a disease context, so stratify or report it separately rather than pooling
> it. This run predates that, so 6579 contributes no histogram, no mode and no scale to any number
> below, and a re-run under the current cohort would not reproduce this shift record.
> 6579 also has no `expression` partition in either root today, so it could not enter regardless
> — `docs/DONOR_6579_REENTRY_2026-09-17.md`.
>
> One further path note: the "retired root" this run reads,
> `core/phenocycler/data/runs/30f6689e4d27d295f890/expression/`, is retired in a stronger sense
> since **2026-09-16 removed the content-addressed run framework entirely**. Stages write to the
> stable `core/phenocycler/data/expression/` and that run tree is historical, though the 20-donor
> cohort output still lives only there.


**The full cohort run completed, but neither new donor has a complete 59-marker alignment.** Donors 6533 and 112 each have 58 markers with finite QC-eligible measurements. `Ker8_18` is acquired and present in both parquet schemas, but its expression is unavailable: every QC-eligible value is missing. The final shift record explicitly marks its scale unavailable for both donors and estimates that marker from the 20 measured donors only.

For each new donor, 46 measured markers pass the declared background-mode SD bar of 0.10. Eleven other measured markers remain above it; these retain their residual flags and cannot support a common cutoff on this criterion. `CD57` remains `refused_reagent_change` and is left unshifted. The outputs are usable for the measured markers subject to these restrictions; they are not ready for a downstream analysis that requires all 59 aligned markers. The cause of the missing upstream `Ker8_18` values is unavailable: this task measured the supplied expression and did not modify the expression pipeline or input files.

**The existing donors change enough to matter.** 15 reference donors changed. Applied scales changed for 296 of 1,157 acquired donor–marker pairs, across 16 markers and all 20 existing donors. The new/old scale ratio ranges from 0.313988504 (`CD11b`, donor 6505) to 1.355924020 (`E_cadherin`, donor 6623). Direct comparison with the saved 20-donor parquet files found 378,780,915 changed cell–marker values. Old absolute thresholds therefore need reassessment before use with this output. The number of changed phenotype assignments is unavailable because downstream typing was not rerun.

**Completed run and inputs.**

The final run started at `2026-09-16T21:16:05.259274+00:00`, finished at `2026-09-16T21:18:21.952646+00:00`, exited with code 0, and took 136.693 seconds. It wrote 22 parquet partitions containing 33,240,038 rows. Histograms used the QC-eligible population (26,665,217 cells across the cohort), finite values at or above 1.0, natural logarithms, and the unchanged 1,024-bin UniFORM registration. The pinned UniFORM commit is `7c33396d12e5a19dc0b36e46a94c3293de0cd814`. The marker set is exactly the 59 markers in the saved 20-donor alignment and Step-1 usability input.

Expression resolves independently for each donor, preferring `core/phenocycler/data/expression/` and warning on fallback to `core/phenocycler/data/runs/30f6689e4d27d295f890/expression/`. The 20 old donors used the retired root; 6533 and 112 used the stable root. `donor_id=0112` resolves to cohort donor `112`, and the output partition and donor column both use `112`. Donors 6457 and 6579 are excluded. The acquisition table records both new donors as Fusion 2.4.1 with 70 channels and panel variant `both`; the resulting `2.4.1/both` group has two donors. The other groups contain 5 (`2.3.1/CD38`), 8 (`2.4.0/CD38`), 4 (`2.4.1/b_Catenin1`), and 3 (`2.4.1/neither`) donors.

| Donor | Input root / partition ID | Output rows | QC-eligible cells | Acquired alignment markers | Markers with finite eligible values |
| --- | --- | --- | --- | --- | --- |
| 115 | retired / 115 | 1,030,803 | 860,863 | 58 | 58 |
| 6374 | retired / 6374 | 2,529,819 | 2,132,181 | 58 | 58 |
| 6380 | retired / 6380 | 737,647 | 623,781 | 58 | 58 |
| 6436 | retired / 6436 | 1,932,797 | 1,579,440 | 58 | 58 |
| 6442 | retired / 6442 | 1,715,283 | 1,419,669 | 58 | 58 |
| 6450 | retired / 6450 | 441,534 | 368,041 | 58 | 58 |
| 6461 | retired / 6461 | 1,877,403 | 1,445,027 | 58 | 58 |
| 6476 | retired / 6476 | 1,979,670 | 1,625,342 | 58 | 58 |
| 6482 | retired / 6482 | 1,366,260 | 1,135,669 | 58 | 58 |
| 6493 | retired / 6493 | 1,545,699 | 1,227,411 | 58 | 58 |
| 6505 | retired / 6505 | 2,909,163 | 2,168,572 | 58 | 58 |
| 6521 | retired / 6521 | 1,399,888 | 1,128,971 | 58 | 58 |
| 6523 | retired / 6523 | 1,023,341 | 811,564 | 57 | 57 |
| 6534 | retired / 6534 | 1,417,265 | 1,145,661 | 58 | 58 |
| 6539 | retired / 6539 | 501,991 | 364,201 | 58 | 58 |
| 6551 | retired / 6551 | 1,845,970 | 1,379,108 | 58 | 58 |
| 6566 | retired / 6566 | 2,208,951 | 1,781,221 | 57 | 57 |
| 6591 | retired / 6591 | 1,367,698 | 1,121,055 | 57 | 57 |
| 6593 | retired / 6593 | 2,046,673 | 1,629,937 | 58 | 58 |
| 6623 | retired / 6623 | 1,576,186 | 1,202,823 | 58 | 58 |
| 6533 | stable / 6533 | 1,027,546 | 879,631 | 59 | 58 |
| 112 | stable / 0112 | 758,451 | 635,049 | 59 | 58 |

The scope of the run is the full donor cohort, with marker availability recorded separately. New-donor `Ker8_18` has 0 finite observations among 879,631 eligible cells for 6533 and among 635,049 for 112; these are measured observation counts, not placeholder intensity values. Its intensity quantiles, donor mode, and shift are **unavailable**, with the reason above. All eligible values are finite for each of the other 58 markers on both new donors. Each acquired alignment marker on every existing donor has measurable histogram support; sparse missing input values are retained rather than imputed.

The existing donors retain 72 missing QC-eligible cell–marker entries across seven donors, in `DAPI`, `FOXP3`, `Ki67`, `PCNA`, `SOX2`, `TCF_1`, `TOX`, and `TP63`. Those individual intensities are unavailable because the input measurement is missing; the remaining observations support each donor's histogram. Per-donor counts are recorded in `verification.json`.

The empty-input guard was necessary: the existing histogram adapter substitutes a one-element log-floor array when a donor has no finite values at or above 1. That produces a synthetic registration target. The runner now excludes only those unmeasurable donor–marker histograms, retains the acquired-donor count, records `unavailable_donors` and an unavailable scale, and preserves the raw missing values in parquet. Marker membership still comes from `panel.donor_markers`; no unacquired marker is added. UniFORM registration itself is unchanged.

**Every reference donor that changed.**

| Marker | 20-donor reference | 22-donor-run reference |
| --- | --- | --- |
| CD107a | 6551 | 6493 |
| CD11b | 6623 | 6534 |
| CD20 | 6551 | 115 |
| CD209 | 6493 | 6523 |
| CD38 | 6482 | 115 |
| CD57 | 6523 | 6380 |
| Caveolin | 6623 | 6533 |
| E_cadherin | 6493 | 6534 |
| HLA_A | 6442 | 6493 |
| M2Gal3 | 6482 | 6533 |
| MPO | 6442 | 6623 |
| PD_L1 | 6521 | 6623 |
| SMA | 6521 | 6623 |
| TP63 | 6623 | 6450 |
| iNOS | 6623 | 6482 |

`CD57` is included in this list because its proposed reference changed, although its shift is refused and every applied scale stays 1.0. Donor 6533 is the new reference for `Caveolin` and `M2Gal3`. Donor 112 is not a reference for any marker.

**Per-marker measurements and the 0.10 bar.**

All numbers below come from the two saved shift records. SD is the population standard deviation of donor histogram-peak positions in natural-log units, using the existing estimator that excludes the first histogram bin. `n` counts donors with measured histograms. `down` means the 22-donor-run SD goes from >0.10 before to ≤0.10 after; `over` means it remains >0.10; `refused` means no shift is applied. Every marker starts above the bar in both runs. Forty-six cross down in each run. **No marker changes its post-registration side of the bar between the 20-donor and 22-donor runs.**

For `Ker8_18`, the 22-donor-run numbers below describe only the 20 measured donors; its complete 22-donor SD is **unavailable**, because the two new donors have no measured modes. For `CD57`, each saved `mode_sd_after` is a diagnostic for the proposed, refused shift. The SD of the actual unshifted output remains 2.346990788 (20 donors) and 2.280491270 (22 donors), respectively.

| Marker | Measured n: 20 → 22 run | Ref 20 | Ref 22 | 20 mode_sd_before | 20 mode_sd_after | 22 mode_sd_before | 22 mode_sd_after | 22 bar result | Provenance in both runs |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| B3TUBB | 20 → 22 | 6534 | 6534 | 0.384905172 | 0.029155610 | 0.367543338 | 0.027926587 | down | automatic_registration |
| Bcl_2 | 20 → 22 | 6623 | 6623 | 0.166536105 | 0.018711168 | 0.163951759 | 0.018586778 | down | automatic_registration |
| Beta_actin | 20 → 22 | 6623 | 6623 | 0.593586842 | 0.339394054 | 0.605607403 | 0.356689291 | over | automatic_registration_over_bar |
| CD107a | 20 → 22 | 6551 | 6493 | 0.267294075 | 0.029452249 | 0.278757214 | 0.029086803 | down | automatic_registration |
| CD11b | 20 → 22 | 6623 | 6534 | 1.097894931 | 0.454250739 | 1.108246319 | 0.627653295 | over | automatic_registration_over_bar |
| CD11c | 20 → 22 | 6482 | 6482 | 0.286193373 | 0.042077182 | 0.415150990 | 0.041785408 | down | automatic_registration |
| CD163 | 20 → 22 | 6482 | 6482 | 0.586481173 | 0.054795063 | 0.560863719 | 0.053222884 | down | automatic_registration |
| CD20 | 20 → 22 | 6551 | 115 | 0.239248588 | 0.035105363 | 0.387497836 | 0.038094353 | down | automatic_registration |
| CD206 | 20 → 22 | 6623 | 6623 | 0.287274439 | 0.017100971 | 0.277718229 | 0.018677024 | down | automatic_registration |
| CD209 | 20 → 22 | 6493 | 6523 | 0.256779854 | 0.016551831 | 0.255087298 | 0.016302510 | down | automatic_registration |
| CD31 | 20 → 22 | 6623 | 6623 | 0.548847854 | 0.154840948 | 0.533225931 | 0.147852684 | over | automatic_registration_over_bar |
| CD34 | 20 → 22 | 6534 | 6534 | 0.353830730 | 0.053561623 | 0.338857355 | 0.055886152 | down | automatic_registration |
| CD38 | 13 → 15 | 6482 | 115 | 0.436315242 | 0.049000031 | 0.455700059 | 0.043813438 | down | automatic_registration |
| CD39 | 20 → 22 | 6493 | 6493 | 0.656550251 | 0.058866810 | 0.650519566 | 0.058586868 | down | automatic_registration |
| CD3e | 20 → 22 | 115 | 115 | 0.312311263 | 0.014452421 | 0.301043680 | 0.015633721 | down | automatic_registration |
| CD4 | 20 → 22 | 6523 | 6523 | 0.198107526 | 0.011990797 | 0.236005500 | 0.011706797 | down | automatic_registration |
| CD44 | 20 → 22 | 6534 | 6534 | 0.709670232 | 0.117788473 | 0.689856317 | 0.123766938 | over | automatic_registration_over_bar |
| CD56 | 20 → 22 | 6493 | 6493 | 0.327455653 | 0.048408275 | 0.327345248 | 0.046831439 | down | automatic_registration |
| CD57 | 20 → 22 | 6523 | 6380 | 2.346990788 | 0.359056640 | 2.280491270 | 0.263023895 | refused | refused_reagent_change |
| CD66 | 20 → 22 | 6623 | 6623 | 0.199847808 | 0.024987199 | 0.348286187 | 0.024762563 | down | automatic_registration |
| CD68 | 20 → 22 | 6591 | 6591 | 0.376367332 | 0.023407855 | 0.367742832 | 0.022484897 | down | automatic_registration |
| CD79a | 20 → 22 | 6482 | 6482 | 0.332108554 | 0.029460568 | 0.349774434 | 0.029055753 | down | automatic_registration |
| CD8 | 20 → 22 | 6523 | 6523 | 0.536521304 | 0.068257558 | 0.796024899 | 0.074588454 | down | automatic_registration |
| CD99 | 20 → 22 | 6534 | 6534 | 0.400673159 | 0.049257541 | 0.453197179 | 0.048105655 | down | automatic_registration |
| Caveolin | 20 → 22 | 6623 | 6533 | 0.508269493 | 0.053553859 | 0.486867300 | 0.099329493 | down | automatic_registration |
| Collagen_IV | 20 → 22 | 6623 | 6623 | 0.204566627 | 0.035097561 | 0.202839145 | 0.034204776 | down | automatic_registration |
| DAPI | 20 → 22 | 6623 | 6623 | 0.279253013 | 0.025197517 | 0.274606056 | 0.024666373 | down | automatic_registration |
| E_cadherin | 20 → 22 | 6493 | 6534 | 0.447494864 | 0.323395331 | 0.547218410 | 0.332844430 | over | automatic_registration_over_bar |
| EpCAM | 20 → 22 | 6623 | 6623 | 0.289425207 | 0.052657543 | 0.311673872 | 0.051198758 | down | automatic_registration |
| FOXP3 | 20 → 22 | 115 | 115 | 0.380154399 | 0.044523465 | 0.365795801 | 0.043537593 | down | automatic_registration |
| GCG | 20 → 22 | 6493 | 6493 | 0.297945553 | 0.067081060 | 0.325067370 | 0.067167813 | down | automatic_registration |
| Granzyme_B | 20 → 22 | 6534 | 6534 | 0.186192430 | 0.020753204 | 0.188307104 | 0.033590522 | down | automatic_registration |
| HLA_A | 20 → 22 | 6442 | 6493 | 0.981445261 | 0.583547037 | 0.988407434 | 0.622866427 | over | automatic_registration_over_bar |
| HLA_DR | 20 → 22 | 115 | 115 | 0.509365104 | 0.029140522 | 0.489419179 | 0.027882145 | down | automatic_registration |
| ICOS | 20 → 22 | 6482 | 6482 | 0.277736399 | 0.019443102 | 0.310808377 | 0.022412727 | down | automatic_registration |
| IDO1 | 20 → 22 | 6482 | 6482 | 0.566498796 | 0.064163612 | 0.549318364 | 0.064360818 | down | automatic_registration |
| INS | 20 → 22 | 6493 | 6493 | 0.717812560 | 0.145394156 | 0.751957126 | 0.149613759 | over | automatic_registration_over_bar |
| Iba1 | 20 → 22 | 6442 | 6442 | 0.308430497 | 0.029051125 | 0.302383800 | 0.027721436 | down | automatic_registration |
| Ker8_18 | 20 → 20 (2 unavailable) | 6534 | 6534 | 0.328257101 | 0.109167623 | 0.328257101 | 0.109167623 | over | automatic_registration_over_bar |
| Keratin_5 | 20 → 22 | 6521 | 6521 | 0.393534054 | 0.035490526 | 0.384691281 | 0.034989986 | down | automatic_registration |
| Ki67 | 20 → 22 | 6523 | 6523 | 0.477518383 | 0.039861838 | 0.469705595 | 0.041109376 | down | automatic_registration |
| LAG3 | 20 → 22 | 6623 | 6623 | 0.381349400 | 0.057356009 | 0.368399480 | 0.056937637 | down | automatic_registration |
| M2Gal3 | 20 → 22 | 6482 | 6533 | 0.556503417 | 0.047038216 | 0.532349722 | 0.044948452 | down | automatic_registration |
| MPO | 20 → 22 | 6442 | 6623 | 0.315447804 | 0.047991089 | 0.310604090 | 0.057549818 | down | automatic_registration |
| PCNA | 20 → 22 | 6493 | 6493 | 0.410303270 | 0.142998344 | 0.425584393 | 0.193486619 | over | automatic_registration_over_bar |
| PD_1 | 20 → 22 | 6623 | 6623 | 0.462658719 | 0.058292995 | 0.444499253 | 0.058899696 | down | automatic_registration |
| PD_L1 | 20 → 22 | 6521 | 6623 | 0.287943576 | 0.023997247 | 0.304676727 | 0.026128877 | down | automatic_registration |
| Pan_Cytokeratin | 20 → 22 | 6534 | 6534 | 0.610029576 | 0.404886167 | 0.588688828 | 0.388870438 | over | automatic_registration_over_bar |
| Podoplanin | 20 → 22 | 6523 | 6523 | 0.467096331 | 0.033012340 | 0.477470340 | 0.036201394 | down | automatic_registration |
| SMA | 20 → 22 | 6521 | 6623 | 0.293246371 | 0.063647314 | 0.324967020 | 0.061417610 | down | automatic_registration |
| SOX2 | 20 → 22 | 6623 | 6623 | 0.455658210 | 0.019736149 | 0.454000930 | 0.019442528 | down | automatic_registration |
| SST | 20 → 22 | 6623 | 6623 | 0.539021224 | 0.067530950 | 0.516393375 | 0.068697912 | down | automatic_registration |
| TCF_1 | 20 → 22 | 6482 | 6482 | 0.382420786 | 0.224733620 | 0.371039125 | 0.214286963 | over | automatic_registration_over_bar |
| TOX | 20 → 22 | 6482 | 6482 | 0.251218229 | 0.017667816 | 0.324537949 | 0.017761963 | down | automatic_registration |
| TP63 | 20 → 22 | 6623 | 6450 | 0.419728537 | 0.015115016 | 0.403774088 | 0.016800481 | down | automatic_registration |
| VISTA | 20 → 22 | 6623 | 6623 | 0.492967725 | 0.099913319 | 0.500717094 | 0.095590603 | down | automatic_registration |
| Vimentin | 20 → 22 | 6493 | 6493 | 1.049877828 | 0.851553297 | 1.050178826 | 0.906152590 | over | automatic_registration_over_bar |
| b_Catenin1 | 4 → 6 | 6623 | 6623 | 0.196047000 | 0.045286162 | 0.227819501 | 0.055212912 | down | automatic_registration |
| iNOS | 20 → 22 | 6623 | 6482 | 0.732057694 | 0.087513642 | 0.785538065 | 0.076800457 | down | automatic_registration |

Median `mode_sd_after` across the 58 automatically registered markers is 0.048199682 in the 20-donor run and 0.047468547 in the 22-donor run. This includes the 20-donor estimate for `Ker8_18` and excludes refused `CD57`. The 12 markers still over the bar are `Beta_actin`, `CD11b`, `CD31`, `CD44`, `E_cadherin`, `HLA_A`, `INS`, `Ker8_18`, `PCNA`, `Pan_Cytokeratin`, `TCF_1`, `Vimentin`. For the new donors, this is eleven measured over-bar markers plus unavailable `Ker8_18`.

**Provenance changes.**

The complete list of markers whose provenance changed is empty: 0 changes.

| Provenance | 20-donor count | 22-donor-run count |
| --- | --- | --- |
| automatic_registration | 46 | 46 |
| automatic_registration_over_bar | 12 | 12 |
| refused_reagent_change | 1 | 1 |
| too_few_donors_to_align | 0 | 0 |

The reagent-change decision uses the unchanged `marker_donor_usability.parquet` measured on the old cohort. `CD57` retains its acquisition-group ratio of 35.121879617 against the refusal bar of 2.0. A fresh 22-donor acquisition-group ratio is **unavailable**, because Step 1 was not recomputed. In particular, `b_Catenin1` has no ratio in the old usability input: its four old donors share one acquisition group, and the estimator requires at least two. Its ratio is recorded as `unavailable` with that input limitation rather than zero. The unchanged provenance therefore describes this run under the retained Step-1 refusal policy, not an independent reassessment of reagent changes in the new panel group.

**CD38 and b_Catenin1.**

`CD38` gains donors 6533 and 112: measured n 13 → 15; reference 6482 → 115. The 20-donor-run mode SD is 0.436315242 → 0.049000031; the 22-donor-run mode SD is 0.455700059 → 0.043813438. Both runs pass 0.10 after alignment and retain `automatic_registration`.

`b_Catenin1` gains donors 6533 and 112: measured n 4 → 6; reference 6623 → 6623. The 20-donor-run mode SD is 0.196047000 → 0.045286162; the 22-donor-run mode SD is 0.227819501 → 0.055212912. Both runs pass 0.10 after alignment and retain `automatic_registration`.

| Marker | New donor | Finite QC-eligible cells | Values ≥1 used in histogram | Applied scale | Aligned log-mode |
| --- | --- | --- | --- | --- | --- |
| CD38 | 6533 | 879,631 | 854,652 | 2.377302350 | 4.298904522 |
| CD38 | 112 | 635,049 | 625,334 | 1.464388240 | 4.298904522 |
| b_Catenin1 | 6533 | 879,631 | 847,006 | 1.349941421 | 7.251478995 |
| b_Catenin1 | 112 | 635,049 | 621,289 | 1.116303120 | 7.381505515 |

The 13 existing CD38 donors have new/old applied scale ratios of 0.865604017–0.892793151. All four existing b_Catenin1 scales are unchanged (ratio 1.000000000). CD38 is absent from output columns for the seven donors that never acquired it; b_Catenin1 is absent for the sixteen that never acquired it. No alignment treats those missing acquisitions as measured zeros.

**Magnitude of changes in the original 20 donors.**

Ratios below are `new applied scale / old applied scale`, restricted to donors that acquired the marker. A changed scale has absolute log-ratio >1e-12, which excludes floating-point noise. These are descriptions of measured changes, not new acceptance thresholds.

| Marker with changed applied scales | Changed / acquired old donors | Minimum ratio | Median ratio | Maximum ratio |
| --- | --- | --- | --- | --- |
| Beta_actin | 19 / 20 | 0.990410184 | 1.000003179 | 1.009698409 |
| CD107a | 20 / 20 | 1.100547930 | 1.111142752 | 1.111142752 |
| CD11b | 19 / 20 | 0.313988504 | 0.866442127 | 1.179965921 |
| CD20 | 6 / 20 | 0.990594465 | 1.000000000 | 1.009494839 |
| CD209 | 20 / 20 | 0.898417016 | 0.898417016 | 0.908092674 |
| CD38 | 13 / 13 | 0.865604017 | 0.874573792 | 0.892793151 |
| CD66 | 19 / 20 | 0.992840733 | 1.001692112 | 1.008862973 |
| Caveolin | 20 / 20 | 0.668003044 | 0.688525546 | 0.817314532 |
| E_cadherin | 20 / 20 | 1.241207654 | 1.303685178 | 1.355924020 |
| HLA_A | 20 / 20 | 0.574624709 | 0.642591090 | 0.732685337 |
| M2Gal3 | 20 / 20 | 1.186932015 | 1.209749309 | 1.244800316 |
| MPO | 20 / 20 | 0.950224279 | 0.959977155 | 0.969830133 |
| PD_L1 | 20 / 20 | 0.848634203 | 0.857970778 | 0.867410072 |
| SMA | 20 / 20 | 0.867875938 | 0.876705094 | 0.885624072 |
| TP63 | 20 / 20 | 1.227687556 | 1.227687556 | 1.241014021 |
| iNOS | 20 / 20 | 0.852179481 | 0.875216307 | 0.889316029 |

Of the 1,157 old donor–marker pairs, 252 differ by more than 1% and 224 by more than 10% in the applied scale. The other 43 markers have identical applied scales for every old donor. Within-donor ranks for each measured marker are preserved by its positive constant scale, but absolute intensities and ratios between differently rescaled markers change.

A changed reference is not the only route to a changed scale: `Beta_actin` retains reference 6623 but the global log minimum changes from 0.00020981262022093053 to 0.0; `CD66` retains reference 6623 but the log maximum changes from 10.756353919289907 to 11.188722915299877. These change the common bin-to-log conversion. Conversely, `CD57` changes its proposed reference but remains unshifted, and six old CD20 donor scales change by roughly one histogram bin despite its changed reference.

**Artifacts, verification, and reproduction.**

The data outputs live in the shared data directory and are not added to Git. The commit contains only the two owned Python modules, the new regression tests, and this report.

- [Shift record](../data/development/celltype/shifts_22donor.json)
- [22 donor parquet partitions and 59 registration figures](../data/development/celltype/uniform22/)
- [Run log](../data/development/celltype/uniform22/run.log)
- [Run manifest with timestamps and source/input hashes](../data/development/celltype/uniform22/run_manifest.json)
- [Input inventory](../data/development/celltype/uniform22/input_inventory.json)
- [New-donor measurement counts and modes](../data/development/celltype/uniform22/new_donor_measurements.json)
- [Full output and baseline verification](../data/development/celltype/uniform22/verification.json)
- [Test results](../data/development/celltype/uniform22/tests.log)

The 22 output parquet files occupy 15,022,621,341 bytes. Verification read all 1,925,108,211 cell–marker entries and checked each output against `input × recorded scale`, or unchanged input where the scale is unavailable, with exact equality including the missing-value masks. Row counts, object IDs and row order, QC flags, coordinates, regions, and islet numbers match the inputs. Donor IDs match the canonical output partition. Baseline parquets also match their original input × baseline scale exactly, with 0 disagreements. Missing values in input rows remain missing; they are not imputed.

The protected `shifts.json` and all 20 existing `uniform/` parquet files (21 files; 14,176,278,038 bytes) have identical SHA-256 hashes, byte sizes, and modification times before and after the run. The acquisition table, Step-1 usability input, and the source files used by the run also match their recorded hashes. The selected tests passed: **25 passed** (19 new tests and six existing UniFORM/baseline tests). The tests cover root preference and fallback warnings, donor aliases, missing and ambiguous inputs, eligibility, marker membership, unavailable measurements, default and custom output isolation, row preservation, and registration of a known synthetic scale factor. Three warnings come from the existing scikit-fda/UniFORM plotting code.

| Artifact | SHA-256 |
| --- | --- |
| Protected 20-donor shifts.json | 6f4e00df17d8d4d507a1393d98a8d4559b33a9e68bea323652e01dfcaee45073 |
| Final shifts_22donor.json | e745cca7d7727f8dbdad5dbcf4cb1546ea494a162d6f45d9f30d27a347bdc92b |
| donor_acquisition_batch.csv | 8cc98a93babc1afbf301c3ec6046021a60289cc2442e116bba784c2d44723c57 |
| marker_donor_usability.parquet | 3282babb6860818108b398095bd853456b71cb0413b9d50e19faf0110131520e |
| panel.py used for the run | 64b803867bf1a4d6167969df4a7cffcac2fa8013ad39c5df4c1e1a65a2fd2213 |
| run_step2.py used for the run | e78e6684a786b4a8e013d1d96e99f14364408d00dde18fc41310f32c7f80bbd1 |

Equivalent execution from the working root is below. This worktree has no local `external/UniFORM` clone, so the invocation supplies the existing shared checkout without changing `uniform_adapter.py`; its commit check still runs. The recorded run additionally writes timestamps and hashes to `run_manifest.json`.

```bash
PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
.venv/bin/python -u - \
  --out data/development/celltype/uniform22 \
  --shifts-out data/development/celltype/shifts_22donor.json <<'PYRUN'
import pathlib, sys
sys.path.insert(0, "core/phenocycler/scripts/celltype")
import run_step2
run_step2.UA.UNIFORM_PATH = pathlib.Path(
    "/home/smith6jt/IO60panc2nd/Islet-Explorer-Senior/external/UniFORM"
)
raise SystemExit(run_step2.main())
PYRUN
```

Omitting `--donors` selects the full eligible cohort. A subset can be named with `--donors`; zero-padded aliases are canonicalized and duplicates removed. Default output paths remain the historical `uniform/`, `figures_step2/`, and `shifts.json`. A custom `--out` isolates parquet, figures, and the default shift JSON; this run explicitly chooses the requested sibling `shifts_22donor.json`. Consumers of this new shift record must honor the `unavailable_donors` reasons and the unavailable scale for new-donor `Ker8_18`; treating it as a numeric scale would be an error.
