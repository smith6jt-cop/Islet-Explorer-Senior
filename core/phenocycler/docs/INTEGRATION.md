# PhenoCycler ↔ Xenium integration

> **Current boundary (verified 2026-09-03).** Geometry is live; identity is
> release-gated. `export_pheno` reads coordinates, areas, QuPath annotation
> context, and identity-safe intensities from immutable core run
> `30f6689e4d27d295f890`. It exported all 20 pancreas sections (25,150,537 cells)
> without inventing a cell type. `import_xenium` likewise ignores embedded
> label-shaped columns unless an independently declared Xenium identity release
> supplies an exact `cell_id` call set. Registration, geometry QC, and geometry
> figures run without identity. Composition, structures, matching, and other
> identity consumers fail before reading caches until the required releases pass
> `phenospace.py`.
>
> The upstream PhenoCycler workflow still stops after expression. Historical
> `phenotype/broad`, `restore_gated_redsea`, and `runs/*/assignments` artifacts
> remain retired and have no compatibility opt-in. Current stage evidence lives
> in `docs/WORKFLOW_2026-09-02_MAP.md` § 4.

Picks up where `python -m phenocycler.pipeline` stops and integrates that with 10x Xenium spatial
transcriptomics processed by the `analysis/xenium` lane (absorbed from `Xenium_Analysis` on
2026-08-28), using the curated panel taxonomy in the `core/panels` subtree (XeniumPanelExplorer).

```bash
# On the workstation the repo .venv is the environment (environment-integration.yml is the HPG recipe).
../../.venv/bin/pip install -e '.[integration]'

python -m phenocycler.integration.manifest --build --report
python -m phenocycler.integration.pipeline --status

# Geometry-only execution is explicit while identity releases are blank.
python -m phenocycler.integration.pipeline --only export_pheno --force
python -m phenocycler.integration.pipeline --only import_xenium --force
python -m phenocycler.integration.pipeline --only register --donor 6436 --donor 6442 --force
python -m phenocycler.integration.pipeline --only qc --force
```

A full unscoped invocation includes identity stages and therefore fails in
preflight while the release paths in `config.ini` are blank. This is deliberate:
selecting a downstream stage or leaving an old output directory in place cannot
bypass provenance validation.

---

## Tissue

The current immutable PhenoCycler run contains **20 pancreas sections**. All 20
pair with Xenium pancreas sections. The Xenium source catalog also contains
pancreatic-lymph-node sections, but this run contains no PhenoCycler pLN section,
so those manifest rows remain `xenium_only`; configured tissue support is not
evidence that a section was measured.

Tissue is a property of the *dataset*, fixed when it is imported, not a mode of the analysis.
It is derived from the ROI (`panc` → pancreas, `pln_1..3` → pancreatic lymph node), written
into the cell table as a column at creation, and validated: a section whose ROI maps to no
known tissue is refused rather than defaulting to pancreas. Every stage therefore runs on one
tissue or on both, and a bare run covers both:

```bash
python -m phenocycler.integration.pipeline                                # both tissues
python -m phenocycler.integration.pipeline --tissue pancreatic_lymph_node # one tissue
python -m phenocycler.integration.pipeline --roi pln_2                    # one section
```

`--roi` beats `--tissue`; with neither, the scope is every ROI of every configured tissue.
Stage runners take that scope and iterate observed manifest rows themselves — they are *not*
called once per ROI, because several output paths are stage-level. Geometry stages run for
sections that exist in the manifest. The PhenoCycler geometry exporter currently refuses any
ROI other than `panc` because the immutable run has no other section to export. Identity
stages preflight the modality-specific release or releases they consume before consulting
cached outputs.

The two tissues diverge in exactly two places, both biological. Everything else — export,
import, vocab, registration, niches, the grid, donor composition, pseudo-cell linking — is
identical.

**1. Structures.** Islets exist in pancreas. A lymph node has none, so hormone-seeded DBSCAN
has nothing to seed on and islet matching (S6) has nothing to match. The structure level is
simply not run there; `--status` marks it `n/a`, and a combined run still matches islets in
the pancreas half. B-cell follicles are the obvious substitute landmark and the code would
support them (`call_landmarks` is generic over a seed lineage), but a follicle called from a
5K RNA panel on one side and a 55-plex protein panel on the other is a far weaker
correspondence than an islet, and shipping it would invite structure-level claims the data
cannot support. `tissues.STRUCTURE_SEEDS` is where it would go.

This propagates into QC. The islet-RMSE gate is *undefined* without islets, not missing, so
for those tissues the verdict rests on tissue Dice alone. Counting the absent metric as a
failed gate would cap every lymph-node donor at `WARN` forever and permanently withhold
`crossmodal` from half the cohort — an artefact of the gate design, not a property of the
registration.

**2. Identity panels.** XeniumPanelExplorer ships `pancreas` and `thymus`, no lymph node. But
the cores a lymph node needs are exactly the *shared* ones — T cell, B/plasma, NK/ILC,
myeloid, granulocyte, endothelial, pericyte, neural, epithelial — and the three
pancreas-specific cores (endocrine, exocrine, fibroblast/stellate) do not apply. So pLN
resolves to the shared pool: 8 cores against the pancreas's 11. That is the correct answer,
not a degraded fallback. Scoring a lymph-node cell against `01_pancreas_endocrine` would not
be harmlessly uninformative — the broad call is an argmax over those axes, so an axis that
cannot exist in the tissue is a competitor that can win.

The same logic narrows donor-level concordance: `Endocrine` and `Exocrine` are dropped from
the comparable-lineage set for pLN, because a near-zero fraction on both sides is two
instruments correctly finding nothing, and letting it into the Spearman would inflate the
concordance with a free correct answer.

`phenocycler/integration/tissues.py` is the single place that knows any of this.

### Section keys on the PhenoCycler side

When both tissues have been acquired, both images are named for the same donor:

```
6539_Scan1.er.qptiff        pancreas
6539pLN_Scan1.er.qptiff     pancreatic lymph node
```

The historical pipeline keyed partitions on `regexp_extract("Image", '[0-9]+')` — the first
digit-run — so both landed in `data/cells/donor_id=6539`. Two sections in one partition is
not a labelling nuisance. Every core stage treats a partition as one image: `redsea.py` reads
`image.iloc[0]` and masks *every* cell in the partition with that one section's GeoJSON,
`restore.py` fits one threshold set per partition, and each section's micron coordinates are
measured from its own slide origin. A shared partition would mask one section with the
other's polygons, normalise two images as one, and interleave two coordinate frames into a
single overlapping point cloud — with nothing raised anywhere.

So the partition key is the **section key**: the donor's digits plus whatever region token
the filename carries, lowercased.

| image | partition | donor | roi |
|---|---|---|---|
| `6539_Scan1.er.qptiff` | `6539` | 6539 | `panc` |
| `6539pLN_Scan1.er.qptiff` | `6539pln` | 6539 | `pln_1` |
| `6539pLN2_Scan1.er.qptiff` | `6539pln2` | 6539 | `pln_2` |

`phenocycler/sections.py` is the single definition, shared by the DuckDB expression that
*builds* the partitions and the Python parser that *takes them apart* — two languages, one
rule, checked against each other in the test suite. An unrecognised token like `6539spleen` is
refused rather than folded into the pancreas, and `redsea.donor_image` asserts its partition
holds exactly one image as a backstop.

A bare `pLN` filename token parses as `pln_1`; numbered tokens preserve their number.
Filename parsing is only a candidate. Export resolves the expression table's exact image id
through `qupath_manifest.json`, then applies exact, reviewable exceptions from the tracked
`data/integration/section_overrides.csv`. That file currently records the HDL115 prefix and
the 6534 `pacn` typo. It may clarify an existing image id, but it cannot add an unmeasured
section or re-key a cell by fuzzy matching.

### The post-Xenium re-stain — protein and RNA on the same cells

After the Xenium run, that section is stained for INS and GCG **on the PhenoCycler**. The
output is an ordinary qptiff, so it flows through the whole core pipeline unchanged. What
differs is the slide underneath — it images the *Xenium* section:

```
PhenoCycler section  ──── serial ────▶  Xenium section
 (55-plex protein)                       (5K RNA)
                                             │
                                        same slide
                                             ▼
                                   post-Xenium re-stain
                                    (INS / GCG protein)
```

This is the only place in the cohort where protein and RNA are measured on **one cell**. The
section key therefore carries a second dimension — which physical slide was imaged — so a
re-stain never lands in the PhenoCycler section's partition:

| image | partition | donor | roi | slide |
|---|---|---|---|---|
| `6539_Scan1.er.qptiff` | `6539` | 6539 | `panc` | pheno |
| `6539xen_Scan1.er.qptiff` | `6539xen` | 6539 | `panc` | **xenium** |
| `6539xenpLN_Scan1.er.qptiff` | `6539xenpln` | 6539 | `pln_1` | **xenium** |

Re-stains export to `cells_pheno_xif/` rather than `cells_pheno/` — same `(donor, roi)`, but
a different piece of tissue, and one partition must stay one section.

`postxen` (S1c) pairs them against the transcriptome. It runs in **both** modes: the
correspondence is a property of that section pair, not of the cohort, so gating it on
`[integration] mode` would answer the wrong question.

> **Measured 2026-08-31 / 09-01 — the "no registration needed" premise is wrong.** The two
> frames sat 31–42 µm apart with a scale mismatch that a median offset cannot remove. The
> re-stain is registered nucleus by nucleus in `scripts/restain/` (0.55–0.63 µm per cell), and
> that lane — not this stage — hosts the protein-vs-RNA measurement; see its README.

**What this unlocks.** INS/GCG/SST are off the 5K panel, so Xenium normally calls Beta/Alpha
from surrogate cores (PDX1, ISL1, NEUROD1, ABCC8). With the re-stain those calls are *scored
against hormone protein on the same cell* — `qc/postxen_qc.csv` reports the agreement and
per-subtype recall. The hormone call is then stamped onto the Xenium table, so S6b types both
sides of the serial pair by protein rather than protein-versus-surrogate-RNA.

**One caveat, enforced by scope.** Post-Xenium tissue has been through protease digestion and
decrosslinking, so epitopes are degraded and intensities are *not* comparable to a fresh
PhenoCycler stain. Only identity travels (`CARRY_COLUMNS`), never raw marker values. CD3e is
excluded on the same grounds — a lower-abundance surface marker does not survive that
chemistry, which is what the reported unreliability predicts.

If a source image uses a one-off noncanonical name, add its exact image id and provenance to
`section_overrides.csv`. A genuinely new naming convention requires a reviewed addition to
`sections.REGION_TO_ROI`; `parse` refuses unknown tokens rather than defaulting and merging
two physical sections.

This also makes pairing exact rather than assumed. The manifest now knows *which* PhenoCycler
sections exist, so a donor whose lymph node was never scanned is `xenium_only` for that ROI
instead of being reported as a pairing that does not exist.

Two consequences worth knowing:

- The core pipeline's `donor_id=*` partitions are **sections**. The current run has 20
  pancreas partitions; a future run containing pLN would add section-keyed partitions rather
  than merge them into those 20. `cfg.discover_sections()` is the work-unit iterator;
  `cfg.discover_donors()` returns unique biological donors.
- `[integration] pheno_tissue_dirs` survives only as an override for a tissue processed
  through a *separate* core-pipeline run. The normal case needs one run and no config.

---

## The two modes, and why the distinction is enforced

```ini
[integration]
mode = sequential      # or same_slide
```

**`sequential`** — serial sections cut from one FFPE block. The two modalities measure
*different cells*. No registration, however good, changes that. Integration therefore
happens at four levels, none of which claims a cell-to-cell correspondence:

| level | pairing unit | needs registration? |
|---|---|---|
| donor | `donor_id` | no |
| niche | kNN composition class | no |
| structure | matched islet | yes |
| pseudo-cell | inferred link | yes, and QC-gated |

**`same_slide`** — one physical section imaged twice. Cells *are* the same cells, so a
genuine paired protein+RNA matrix is recoverable and `sameslide.py` produces it.

`pair_donor` takes a `method`:

| `method` | criterion |
|---|---|
| `iou` | **polygon overlap** — two segmentations of one cell overlap substantially, two neighbours do not |
| `centroid` | mutual-nearest centroid within a 5 µm cap |
| `auto` (default) | IoU, falling back to centroids when polygons are unavailable |

IoU is preferred wherever it can run. In packed islet tissue an adjacent nucleus is often
nearer than the partner outline's centroid, which is precisely the case centroid distance
gets wrong and overlap does not.

Polygons come from `boundaries.py`, which reads QuPath GeoJSON (`data/redsea_scratch/geojson/`)
and Xenium `cell_boundaries.parquet`, and warps the moving side through the transform
**vertex-wise** so the non-rigid field is applied, not just its affine part.

The two formats disagree about units — the GeoJSON is in full-resolution qptiff **pixels**,
every contract table is in **microns** — and getting that wrong does not error, it puts one
section's polygons a thousand microns from their own cells. Rather than read the qptiff's
resolution tags (which needs the image mounted and returns nothing when the tags are absent,
so a missing scale silently becomes 1.0), the conversion is **fitted from the data**: both
files describe the same cells under the same ids, so a per-axis linear fit of polygon centroid
against `x_um`/`y_um` recovers it exactly. A fitted scale outside 0.05–2.0 µm/px, or a median
residual above 5 µm, means the two files are not the same section and is an error rather than
a silently wrong overlay.

`auto` records which matcher actually ran in `stats['method']`. A silent fallback would make a
polygon-loading failure look like a successful centroid run, and the two have different error
characteristics. `method='iou'` re-raises instead of falling back, so an explicit request
cannot be quietly downgraded.

The same matching is used by `postxen.py` (S1c) for the post-Xenium re-stain, where it matters
most: that really is one section imaged twice, so overlap is the correct criterion by
construction.

Every mode-specific module refuses to run in the wrong mode, as a hard error rather than a
warning. The validity of every downstream claim rests on which situation holds, and a
warning is too easy to scroll past.

---

## Pipeline

```
manifest      S0   donor ↔ Xenium-run pairing; stamps donor_id/roi into the data
export_pheno  S1a  immutable-run geometry → cells_pheno/; optional exact release attachment
import_xenium S1b  h5ad | SpatialData zarr | raw bundle → cells_xen/; labels release-only
postxen       S1c  post-Xenium re-stain pairing                       [Xenium identity]
vocab         S2   release-declared lineage vocabulary + marker/gene crosswalk
structures    S3   islets / ducts / vessels                           [identity, pancreas]
register      S4   coarse rigid → affine → optional non-rigid         [geometry]
match         S6   islet ↔ islet assignment                           [identity, sequential]
cellmatch     S6b  endocrine cells within matched islets              [identity, sequential]
sameslide     S11  cell ↔ cell assignment                             [identity, same_slide]
grid          S7   identity-composition niches + registered hex grid  [identity]
donor         S8   donor-level composition concordance                [identity]
crossmodal    S9   pseudo-cell linking                                [identity, inference]
qc            S10  registration gates; omit identity sections when unavailable
figures       S10  geometry overlays; add composition only with releases
```

Everything after S1 speaks one modality-agnostic parquet schema (`contract.py`). Parquet
metadata stamps schema version 1, the `geometry` or `identity` profile, and the identity
release id (blank for geometry). Nullable identity columns in a geometry table mean
unavailable; they are never coerced to negative evidence. An identity table is produced only
by exact native-key attachment (`object_id` UUID for PhenoCycler, `cell_id` for Xenium), with
equal, duplicate-free cell sets and a release-declared label space.

---

## What the data actually looks like

### The join key had to be built, not just read

`donor_id` is a clean key on paper: both repos already share `donor_metadata_panc.xlsx`
verbatim, with the same columns and the same nPOD donors (26 appear in the vendored
`xenium_paths.csv`; the run cohort is a 20-donor subset of them). But:

- **No Xenium artifact carries it.** Every Xenium output is keyed by the XETG slide serial
  (`0041323`), which encodes no donor, no region and no disease status.
- **`xenium_paths.csv` — the one file that maps donors to runs — is read by nothing** in the
  Xenium repo. It is hand-maintained. It is vendored here at
  `data/integration/xenium_paths.csv` and made load-bearing.

  It is no longer *unenforced*: `--verify-paths` resolves every recorded bundle path on the
  current filesystem, applying `xenium_root`, and exits non-zero naming each one that is
  missing, is not a bundle, or was left unrewritten. Because donor identity exists nowhere
  else — the bundles are keyed by serial — this file cannot be regenerated from the data.
  Repair it; do not replace it.

  As vendored (2026-08-14) it holds **56 rows covering 26 donors**, including all 20 of the
  current run cohort — which is what settles "which 20 of the 26". 14 all-blank spreadsheet
  spacer rows were stripped on that date. One genuine gap remains: **donor 6591 `pln_1` has
  no recorded path**, and is now reported rather than silently dropped.
- **The processed samples are not in it.** `0041323` and `0041326` have finished h5ads and
  appear in no manifest, so their donor is recorded nowhere. They enter the manifest as
  `pair_status = donor_unknown` and are excluded from paired analyses until someone fills in
  `data/integration/donor_overrides.csv`. This is visible in `--report`, not silent.
- **A slide serial is not a section.** One serial can carry both a pancreas and a
  pancreatic-lymph-node region (donor 6539 → `0059865` with `__Panc__` *and* `__pLN__`
  bundles), so `data/raw/{serial}.zarr` collides. The section key here is
  `{serial}__{roi}`.

`donor_id` is normalised to a string everywhere — the workbook stores it numerically,
`xenium_paths.csv` as text, PhenoCycler derives it by regex. Left alone, `6539` and `6539.0`
are different keys and the merge yields zero pairs.

### Identity vocabularies are release declarations

```
PhenoCycler release   producer labels ── declared label_space.json ──▶ common lineage
Xenium release        producer labels ── declared label_space.json ──▶ common lineage
```

There is no module-global PhenoCycler phenotype roster. Each release declares every producer
label, whether it maps at the common resolution, inherits from a declared parent, or is
explicitly `Ambiguous`/`Unavailable`. Undeclared labels are errors. The former broad mapping
exists only in `integration/data/legacy_broad_space.json` for explicit audit tests; loading it
does not make it a current identity source.

The Xenium panel/gene half of the crosswalk is *derived* from the `core/panels` subtree
(XeniumPanelExplorer @ `076d1192`,
verified byte-identical 2026-08-28), which ships `panel_roles.csv` and
`identity_core/*.csv` specifically as a machine-readable contract for Python consumers. A
drift guard fails loudly if an upstream lineage stops being covered.

### The `evidence_positive` correction

`evidence_positive` is nullable in geometry tables and required to be fully observed before
composition or concordance is computed. An identity release may mark producer fallbacks as
false, but unavailable evidence is not a negative call. Every composition figure reports
`frac_all`, `frac_evidence_positive`, and `frac_fallback`; none is emitted without declared
identity provenance.

### Off-panel markers

The deployed 5,081-gene panel omits **INS, GCG, SST and Vimentin** — markers historically
used for PhenoCycler endocrine and fibroblast calls. Any declared identity comparison must
therefore record Xenium's curated surrogate identity panels (PDX1, ISL1, NEUROD1, ABCC8 …
for endocrine) rather than pretending those are 1:1 measurements.

That is orthogonal evidence rather than circular confirmation — scientifically a feature —
but it thins the anchor set available for pseudo-cell linking, so `crossmodal.py` refuses to
run below `crossmodal_min_anchors` rather than producing links from three genes.

---

## Registration

Two candidate tracks are available; every track with the required inputs runs and the
best-by-QC choice is recorded.

**Track A — point set.** RANSAC over declared structure correspondences → ICP over mutual
nearest neighbours. It needs identity-derived structures, so it is unavailable for the
current geometry-only cohort even though it needs no images.

**Track B — raster/image.** Coarse orientation search → opencv ECC (Euclidean then affine) →
optional SimpleITK B-spline. Uses DAPI↔DAPI where available and an all-cell density raster
otherwise; it is the current geometry-only path.

Three things that matter:

- **Rotation and mirroring are not optional.** Sections are mounted by hand at arbitrary
  rotation and are sometimes flipped. The coarse stage searches rotation exhaustively and
  evaluates the mirror hypothesis explicitly; `fit_rigid_umeyama(allow_reflection=True)` can
  represent the flip, which textbook Umeyama cannot.
- **Non-rigid deformation is real but bounded.** The B-spline is coarse (control points
  hundreds of microns apart) and hard-capped at `reg_max_disp_um`. Given enough freedom a
  B-spline will align *any* two images, and a "successful" registration of unrelated
  sections is worse than an obvious failure.
- **Resolution is managed explicitly.** The coarse search runs at ~25 µm/px, ECC and the
  B-spline are capped by pixel count, and QC runs at 10 µm. Because `Transform` holds its
  affine in **microns**, a transform estimated at any resolution drops into the next stage
  unchanged. Registering a 10–20 mm section at the 2 µm working resolution would otherwise
  be intractable.

**Reading morphology images.** The Xenium ingest sets `morphology_focus=False`,
`morphology_mip=False`, `aligned_images=False` ("skip the 35GB TIFFs"), so the zarr does not
contain them. `rasterize.read_xenium_morphology()` reads `morphology_focus/` directly from
the original bundle, at a pyramid level near the target resolution.

### QC gates

| verdict | condition | permits |
|---|---|---|
| PASS | tissue Dice ≥ 0.80 and, when identity structures exist, islet RMSE ≤ 150 µm | everything allowed by available profiles |
| WARN | one gate met | all but pseudo-cell linking |
| FAIL | neither | donor-level and niche only |

When the tissue has no structural unit, or when identity is not attached and therefore no
valid structure can be called, the islet-RMSE gate does not exist. The verdict is PASS on
Dice alone and FAIL otherwise. The dropped gate and reason are written explicitly.

Verdicts are per **(donor, ROI)**, not per donor: a donor's pancreas can register cleanly
while its lymph node does not, and `crossmodal` gates on the section it is about to link.

A FAIL donor is **not** dropped from the study — only from the spatial analyses. Losing good
donor-level data because two sections would not align would be the wrong trade.

---

## Outputs

```
data/integration/
├── manifest.csv                     pairing table + pair_status + tissue [tracked]
├── vocab_crosswalk.csv              lineage + marker crosswalk, all tissues [tracked]
├── vocab_crosswalk_<tissue>.csv     per-tissue crosswalk (cores differ)  [tracked]
├── xenium_paths.csv                 vendored donor → run map           [tracked]
├── donor_overrides.csv              donor ids not recorded upstream    [tracked]
├── section_overrides.csv            exact PhenoCycler image exceptions [tracked]
├── cells_pheno/donor_id=*/roi=*/    schema/profile-stamped contract tables
├── cells_xen/donor_id=*/roi=*/
├── structures/{phenocycler,xenium}/donor_id=*/roi=*/{islets,ducts,vessels}.parquet
├── registration/donor_id=*/roi=*/   transform.json, displacement.npz, qc.json
├── paired/islets.parquet            matched islet pairs      [sequential]
├── paired/cells.parquet             paired protein+RNA cells [same_slide]
├── niches/                          niche_profiles.csv, grid_bins.parquet
├── crossmodal/pseudocell_links.parquet   INFERRED links
├── qc/                              registration_qc.csv, composition.csv,
│                                    disease_trend.csv, integration_report.txt
├── figures/                         registration overlays, concordance, trends
└── export/                          QuPath + Xenium Explorer round-trips
```

The Xenium Explorer export uses the **`cell_id`-matched CSV**, never the zarr path: the zarr
writer is position-indexed and the Xenium pipeline drops ~7% of cells at QC, so positions no
longer align with `cells.parquet`.

---

## Verifying it works

```bash
pip install -e '.[integration,test]'
pytest -q tests/integration                          # explicit optional lane
python -m phenocycler.integration.pipeline --status
python -m phenocycler.integration.manifest --build --report

# Do this FIRST on any machine that actually mounts the Xenium tree. Everything
# downstream is worthless if the manifest points at paths that do not resolve.
python -m phenocycler.integration.manifest --verify-paths
python -m phenocycler.integration.manifest --verify-paths --donor 6436 --donor 6442
python -m phenocycler.integration.pipeline --only export_pheno --donor 6436 --donor 6442 --force
python -m phenocycler.integration.pipeline --only import_xenium --donor 6436 --donor 6442 --force
python -m phenocycler.integration.pipeline --only register --donor 6436 --donor 6442 --force
python -m phenocycler.integration.pipeline --only qc --force

# Expected to halt before reading caches until both identity releases are configured:
python -m phenocycler.integration.pipeline --only structures --donor 6436
```

The test suite asserts recovery of *known* transforms — a synthetic warp is applied and the
registration must recover it, including the mirrored case and the non-rigid displacement
sign. Real affine results retained for 6436 and 6442 have tissue Dice 0.586 and 0.768,
respectively, and both correctly fail the fixed 0.80 gate. A 2-million-pixel, 200-iteration
real B-spline probe on 6436 was stopped after extended non-completion; it is performance
evidence, not a completed registration result.

**Scientific acceptance**, in order of how much they demand:

1. Donor composition should correlate across modalities (evidence-positive fractions).
2. The ND → AAB → T1D β-loss gradient should appear in **both** modalities independently.
3. On matched islets, protein INS⁺ fraction should track Xenium Beta fraction, and insulitis
   grades should agree.
4. Permuted-pairing and rotation nulls must destroy all of the above.

## Notes

- The core pipeline is untouched: `python -m phenocycler.pipeline` still runs its four
  live stages (ingest → geometry → redsea → expression) with no new dependencies.
- `Xenium_Analysis/notebooks/06_xenium_phenocycler_integration.ipynb` is superseded. It
  applied a 40×/20× "magnification" scale factor to coordinates that are already in microns,
  synthesised a disc segmentation instead of using the real polygons in the zarr, registered
  with translation only, and was never wired to the real cohort.
- `Xenium_Analysis` was absorbed as `analysis/xenium` on 2026-08-28 with its 39 hardcoded
  `/blue/maigan` paths cleared onto two roots (`ROOT` for code, `XENIUM_DATA_ROOT` for data).
  This package still reads its *artifacts* rather than importing it; the handful of shared
  constants (DBSCAN parameters, insulitis thresholds) remain vendored with provenance comments.
