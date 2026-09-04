# Historical handoff — region layer complete, 2026-08-11

> Consolidation note (2026-09-03): this preserves the source-repository
> analysis record. Runtime paths and identity rules are superseded by
> `../README.md`, `../CLAUDE.md`, `../paths.py`, and
> `../config/sample_manifest.csv`. Never recover donor identity by parsing a
> filename, and never assume large data live below the tracked code tree.

State of the SpleenFollicleCounterQP analysis after the region layer was
rebuilt, exported in full from QuPath, and wired into every notebook.

Supersedes the handoff of earlier the same day, which recorded the QuPath work
as blocked. It was not — see below.

---

## Where the project stands

The **region layer** — how every object is assigned to Follicle / PALS /
Perifollicular / RedPulp / Trabeculae — is the single geometry the whole
analysis reads. It is built from the QuPath pixel-classifier exports, and
nothing resolves regions through the annotation hierarchy or through QuPath's
signed-distance columns any more.

**Headline result: rs3184504 is null on every region.** SmallVessel density,
per donor, n=14:

| Region | Kruskal-Wallis | Spearman dosage |
|---|---|---|
| Follicle | p = 0.724 | ρ = +0.081, p = 0.782 |
| PALS | p = 0.318 | ρ = −0.413, p = 0.142 |
| Perifollicular | p = 0.519 | ρ = −0.294, p = 0.308 |
| RedPulp | p = 0.772 | ρ = −0.069, p = 0.814 |
| Trabeculae | p = 0.997 | ρ = −0.012, p = 0.968 |

At 4 / 7 / 3 donors only very large effects are detectable, so this is weak
evidence of absence rather than evidence of no effect.

These moved against the earlier handoff, which reported Follicle ρ = −0.017 and
had no Trabeculae row at all. Three reasons, all corrections: 1901HBMP004's
refreshed geometry took it from 450 to 811 Follicle-resident vessels, which
flipped the sign of the Follicle trend; Trabeculae now exists for all 14 donors
rather than 3, so it is reportable; and its area left RedPulp. PALS and RedPulp
are unchanged.

**`scripts/compute_primary_endpoint.py` writes these tables.** Nothing did
before — `plot_region_layer_results.py` only *reads* them, so after the layer
changed underneath, the figures kept drawing the old numbers and nothing failed.
Run it before the plotting script, always.

### One nominally-significant result, and why it is an artifact

Moving H4 to n=14 produced the project's only p < 0.05 dosage trend:
**Follicle/PALS area ratio, Spearman ρ = +0.561, p = 0.037**. Do not report it.

Two donors carry it entirely. Every donor sits between 0.30 and 2.25 except
`HDL073` at **10.1** and `1901HBMP004` at **19.8** — both T/T, both CODEX. The
ratio explodes because its denominator nearly vanishes: PALS is 1.6% and 1.1% of
tissue in those two, against 3–16% elsewhere, from 74 and **5** annotated
objects against 61–2233 in the Phenocycler donors.

| stratum | n | ρ | p |
|---|---|---|---|
| all | 14 | +0.561 | **0.037** |
| Phenocycler only | 9 | +0.224 | 0.563 |
| CODEX only | 5 | +0.949 | 0.014 |

Within Phenocycler there is nothing. The signal lives entirely in the CODEX
stratum, where n=5 and the same two donors dominate. A Mann-Whitney on platform
misses it (p = 0.797) because the ratio is so skewed that medians hide the two
outliers — so the standard platform check does **not** catch this one.

It is also not significant after correction: p = 0.037 is one of 8 metrics in
H4, giving BH q ≈ 0.29.

The sparse PALS is upstream, not a region-layer defect — PALS reproduces the
annotation CSV at 0.93–1.00× for every donor, including these two. The
classifier simply annotates very little PALS on the CODEX sections. That makes
any *ratio* with PALS in the denominator fragile by construction, which is worth
remembering before adding another one.

**The design is identifiable in the annotation cohort.** CODEX 1/2/2 vs
Phenocycler 3/5/1 across C/C, C/T, T/T gives corr(platform, dosage) = +0.288,
VIF 1.09, design matrix full rank (`region_layer_identifiability.csv`).

---

## What changed in this round

### The QuPath work was never blocked

QuPath **0.7.0** — the version `project.qpproj` is saved with — is installed at
`/blue/maigan/ebutterworth/software/QuPath/bin/QuPath`, is group-executable, and
runs headless with no display. The previous session checked only `$HOME/QuPath*`
and `which qupath`. There is also a `module load qupath/0.4.3`, too old to open
a 0.7.0 project.

```bash
QP=/blue/maigan/ebutterworth/software/QuPath/bin/QuPath
"$QP" script -p "$PWD/project.qpproj" "$PWD/scripts/export_all_regions_geojson.groovy"
```

Give `--project` an **absolute** path. QuPath derives the project's data
directory from the file's parent, so a bare `project.qpproj` dies with
`NullPointerException: … "this.dirBase" is null`.

`DnnTools` fails to load (`GLIBCXX_3.4.30 not found`). That disables OpenCV DNN
— pixel classifiers and InstanSeg — so the segmentation itself cannot be re-run
here. The annotation export needs neither.

### Trabeculae, and the 1901HBMP004 flag

Every donor in the project had trabecular annotations all along; the old
exporter filtered to Follicle + PALS. `trabeculae_source` is now `qupath` for
**14 of 14** donors, not 3, and that area has moved out of RedPulp (HDL063
182.3 → 173.4 mm²).

`1901HBMP004`'s export was stale — 54 Follicle and 28 PALS objects against the
35 and 5 in the CSV, reproducing 0.82× / 1.59× of the CSV areas. Re-exported it
matches exactly and reproduces 0.97× / 0.99×. **Only HDL073 is still flagged**,
at 0.86×, which is real: it has many small follicles, so the fragment filter
removes more of its area than the 2.6% cohort average.

### Region assignment is point-in-polygon everywhere

`Parent` did not fail randomly. It failed for one donor, stranding 21.5% of
`1901HBMP004`'s SmallVessels at the image root against a cohort norm under 0.3%
— and that donor is one of only three T/T. A per-image failure of that size in a
14-donor design is a genotype effect in all but name. It reported **0**
Follicle-resident vessels there; containment finds **811**.

- `get_vessels`, `compute_density`, `assign_vessels_to_follicles` all take the
  polygon route. `get_vessels(..., from_parent=True)` reproduces old numbers.
- `compute_density` delegates to `compute_density_from_regions`, so areas come
  from polygon unions. Summing annotation areas double-counted the CODEX
  `Red_Pulp`/`Sinusoid` pair.
- H11/H15/H16/H17/H20 assign by containment instead of taking an argmin over
  QuPath's signed-distance columns, which were measured against the pre-rebuild
  geometry.
- `signed_distance_by_image` measures Follicle and PALS against the same
  polygons the labels came from. QuPath's columns measure to the *unfiltered*
  annotations, so a cell inside a removed fragment read as inside a follicle
  while being labelled RedPulp.
- **SmallVessel/LargeVessel/Trabeculae distances stay with QuPath.** The region
  layer holds no vessel polygons by design, so H17's vessel rings are unchanged;
  only its tissue-compartment labels moved.

### Reproducibility gaps closed

- `radial_pfz_width.csv` had **no producer** despite the docstring claiming one.
  `radial_follicle_profiles.py` now derives it — SMA peak, then the half-max
  crossing above the red-pulp baseline — and reproduces the on-disk values
  exactly for all 14 donors.
- `fit_area_floors.py` derives the GMM area floors, which had no fitting code.
- `scripts/build_region_layer.py` is deleted. It was the abandoned cell-derived
  builder and defaulted to the **same** `--out-dir` and `--table` as the current
  one, so running it silently overwrote the live layer with a different schema.
- `tests/test_region_layer.py` — 38 tests, no large data needed. Full suite 117.

**The PALS floor is not what its comment claimed.** A 2-component GMM on PALS
log-area lands at 4,141 µm² and would discard 67% of objects and 9.5% of PALS
area, contradicting the validation that requires PALS to reproduce the CSV at
0.95–1.00×. 612 µm² is a geometric floor — roughly one arteriole's
cross-section. The Follicle floor does fit: 8,490 against the committed 7,945,
close enough that neither the object count nor the area dropped moves
materially. Both constants stand; the comments now say why.

### The annotation notebooks were on two different cohorts

H2, H3, H4, H5 and H6 called `load_data()` — AnnotationsFinal.csv alone, **11
donors** — while H1, H7, H8, H9, H10 and H11 called `load_all_data()`, which
adds the three ForSH2B3 CODEX donors for **14**. So follicle structure, tissue
proportions, vessel morphology and follicle vascularization were computed on a
different, smaller cohort than vessel density, and nothing recorded it.

The two CSVs carry identical column sets, so nothing forced the split; H1/H7
were migrated at some point and the other five were missed. All now use
`load_all_data()`. Their numbers move — three donors added, two of them C/T.

**The n=11 cohort had the same collinearity the single-cell cohort has.** Drop
the three ForSH2B3 donors and the only CODEX donors left are `1901HBMP004` and
`HDL073` — *both T/T*. So platform was collinear with the T/T contrast in H2–H6
exactly as it is in H14–H17/H20, and no platform-adjusted reading of those five
was identifiable. Going to n=14 fixes it: CODEX 1/2/2 vs Phenocycler 3/5/1. It
is the same defect as open item 1, and worth checking for anywhere a notebook
reads a narrower source than the cohort definition.

Worth checking the same way if a notebook ever looks out of step: `load_data()`
is not the cohort, `load_all_data()` is.

### No analysis crops any more

H12, H13, H14, H15, H16 and H17 cropped the Phenocycler sections to a
CODEX-sized window (~5000 µm, or 5000×3000) before analysing them. H20 had
already stopped; the rest had not followed. All six now analyse full sections.

The reasoning was already in `PLATFORM_CONFOUND.md` §5 and simply had not been
acted on: cropping discards 75–88% of the data to address only the first of the
four confound mechanisms — unequal effort — and leaves field selection, optical
resolution and detection sensitivity untouched. A larger section gives a more
precise rate, not a biased one.

The thresholds are set to `np.inf` rather than deleted, so the window machinery
still runs and returns "no crop" for every sample; the downstream guards stay
exercised instead of becoming dead branches. Verified by replicating both window
rules against `AnnotationsFinal.csv` — no sample is cropped under either, and
the largest extent in the cohort is 23,490 µm.

**This leaves H13 misnamed.** Its whole identity was being the "cropped
replication" of H1/H2/H5. What still distinguishes it is the deliberate
re-inclusion of HDL018/HDL021, so it is now an excluded-donor sensitivity
analysis under the wrong title. Rename or retire it.

### analysis/figures/ is wiped before a full run

The directory is gitignored, so stale figures accumulate invisibly and nothing
in git shows it. On 2026-08-11 it held 95 PNGs spanning several partial runs.
Delete `analysis/figures/*.png` before a full re-run, then confirm afterwards
that every file is newer than the wipe — otherwise a figure whose producing cell
was removed simply stays there looking current.

### numpy 2.0

`ndarray.ptp()` is removed in 2.0. 14 call sites across H12–H17 moved to
`np.ptp()`.

---

## Verified state at the end of this round

`scripts/verify_analysis_state.py --since '<run start>'` is green on every
check. It exists because each thing it tests failed silently at least once here.

- **20 of 20 notebooks execute cleanly.** 139 figures, 126 tables, all current.
- **118 tests pass.**
- No notebook resolves a region through `Parent`, takes an argmin over QuPath
  signed distances, reads the superseded `analysis/geojson/`, carries a local
  pixel-size table, calls the removed `ndarray.ptp()`, or can crop.
- Region layer: 14 donors, Trabeculae for all of them, only HDL073 flagged.

The cell-derived boundary method is now quantified rather than just asserted to
have failed. H14 compares it against the region layer: mean Follicle IoU
**0.53**, mean PALS IoU **0.18**. Two markers do not reproduce a
classifier trained on 25-33 channels, and PALS is where it collapses.

## Open items

### 1. Take the single-cell cohort to n=14 — the highest-value item

H14/H15/H16/H17/H20 read `Cells.csv` directly, so they cover **11 donors**, in
which *both* CODEX donors are T/T and platform-adjusted genotype effects are
**not identifiable**.

That is now an ingest artefact, not a data limit. `HBMP006`, `HBMP007` and
`HDL098` *do* have single-cell data — 289,500 / 284,836 / 205,792 cells in
`Measurements/Spleen_Upload2_2024Sept17_v2_distances.csv`, which
`load_single_cells()` already merges and the region-layer build already
consumes. What kept them out was that this CSV carries
`Signed distance to annotation Sinusoid/Vessel µm` and **not** the
Follicle/PALS/RedPulp columns those notebooks took an argmin over. Region
assignment no longer needs them.

Moving those notebooks from `pd.read_csv(CELLS_CSV, …)` to `load_single_cells()`
would give the single-cell arm the same identifiable design the annotation
cohort has. The work is marker-name harmonisation, which `load_single_cells`
already does, plus checking which of the 22 clustering markers survive in the
new sections (`MARKERS_LOST_IN_NEW_CODEX` lists 7 that do not).

### 2. Environment traps

- This shell exports **`MPLBACKEND=Agg`**. Under Agg, `plt.show()` silently
  emits nothing inside a Jupyter kernel: figures still reach
  `analysis/figures/` through `save_figure`, but none is embedded in the
  `.ipynb`. Nothing errors. `scripts/run_notebooks.sh` forces the inline
  backend; do the same for any ad-hoc run.
- `pytest` had to be installed into `xenium_analysis`.

### 3. Known data defects, unchanged

- **Markers are background-subtracted with negatives clipped at zero.** No
  marker anywhere has a negative value; the exactly-zero fraction runs 0.1% in
  Phenocycler donors to **86%** for CD3e in 20_06_CC2B. Any fixed-percentile
  marker gate is platform-dependent by construction — do not reintroduce one.
- **Only 18 of 39 markers have data in every donor**
  (`MARKERS_USABLE_ALL_DONORS`). CD35 and Lyve1 are entirely NaN for all 9
  Phenocycler donors; CD11c is all-NaN in HBMP007. Check per-donor non-null
  fractions, never column presence.
- `1901HBMP004`'s annotation hierarchy is still broken in the CSV. Nothing reads
  it now, so it is no longer worth a QuPath re-run.
- 4–9% of SmallVessels fall outside every region polygon, just beyond the KDE
  tissue contour at section edges. They are labelled `"Unassigned"` rather than
  dropped. H1 warns above 10%.

---

## Environment

```bash
PY=/blue/maigan/smith6jt/miniforge3/envs/xenium_analysis/bin/python   # 3.12.13
QP=/blue/maigan/ebutterworth/software/QuPath/bin/QuPath              # 0.7.0
```

Run notebooks with cwd = `analysis/`. `xenium_analysis` is the **only** env with
full coverage — KINTSUGI lacks scanpy/anndata/harmonypy/leidenalg and cannot run
H15/H16/H17, despite older plans naming it.

Rebuilding the layer from nothing (all of `analysis/regions/` is gitignored):

```bash
"$QP" script -p "$PWD/project.qpproj" "$PWD/scripts/export_all_regions_geojson.groovy"
$PY scripts/radial_follicle_profiles.py      # -> radial_pfz_width.csv
$PY scripts/build_region_layer_qupath.py     # -> analysis/regions/ + QC table
$PY scripts/validate_region_layer.py
$PY scripts/compute_primary_endpoint.py     # -> P1_vessel_density_*.csv
$PY scripts/plot_region_layer_results.py    # reads those, draws P1_*/P2_*
$PY scripts/build_region_layer_qupath.py --include-excluded \
      --images HDL018 HDL021 --table analysis/tables/region_layer_qc_excluded.csv
bash scripts/run_notebooks.sh
```

The last build is for H12/H13 only, which deliberately re-include HDL018 and
HDL021 to enlarge the C/T group. HDL018 reproduces its Follicle area at 0.83× —
worth remembering when reading those two.

Pixel sizes are **not** readable from the OME-XML. Use `get_pixel_size`:
HBMP006/HBMP007/HDL098 = 0.37740007578193524; 1901HBMP004 and HDL073 = 0.374
(their headers wrongly declare 0.5078, a stale value written by
`convert_to_ome_tiff.py`); everything else 0.5077663810243286.

---

## Approaches that were tried and failed

Recorded so they are not repeated.

- **Deriving Follicle/PALS from CD20/CD3e alone.** Two markers cannot match a
  classifier trained on 25–33 channels. Comparing z-scored marker fields labels
  anything below CD20's mean as PALS whether or not CD3e is elevated — HDL098
  got a "PALS" whose CD3e median (40) was *lower* than its own red pulp (71).
- **Rank/ECDF marker normalization.** Discarding magnitude makes spatial contrast
  a function of the clipped-zero fraction, so the marker with more clipping wins
  the dominance test. It inverted HDL073, a visibly follicle-dominated section,
  to PALS-dominant.
- **Otsu thresholds**, on either intensity or the density field. On intensity it
  calls 40–73% of cells positive; on the combined field it swings 1.22–2.87
  between donors and at the high end erases the follicles outright.
- **A QC gate that zeroes an oversized region.** This is what emptied PALS for
  HBMP007 and HDL098. An empty PALS is the method reporting its own failure —
  every spleen has a T zone. Flag, never drop.
- **Anchoring the perifollicular band outside the PALS.** Wrong by definition;
  and anchoring it on all white pulp instead traces the tissue edge, because
  PALS is thousands of objects at median 55 µm across.
- **A GMM crossover as the PALS area floor.** See above — it deletes 9.5% of
  PALS area and breaks the CSV agreement.

---

## Provenance

The pre-rebuild state is archived at
`/blue/maigan/smith6jt/SpleenFollicleCounterQP_archive/2026-08-05_pre-corrected-data/`
— 256 files, checksummed, tagged `pre-corrected-data-2026-08-05` (commit
`37f1cf4`).

**Absolute numbers changed scale during this rebuild.** Do not compare
post-rebuild follicle areas, counts or densities against the archived ones.
