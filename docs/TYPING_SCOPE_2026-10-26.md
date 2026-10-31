# Stage 06 typing scope: both assays and same-slide acquisitions

## Owner typing updates, 2026-10-27

The CODEX Islet Explorer's committed work (`205be40`, `2d0bce1`, `d032f29`) and the
maintainer's later tie clarification also apply to the Senior protein path:

- If a nucleus fills a segmented cell and REDSEA's selected cytoplasm value is
  undefined, read that marker's measured whole-cell value. Preserve a measured
  zero; if the whole-cell value is also missing, leave the measurement missing.
  A cell without a usable nucleus remains excluded from analysis.
- Establish the broad T-cell lineage first. Keep any specific rule call. For a
  T cell still lacking a subtype, compare CD8 and CD4 in the same normalized
  units within its donor, before clipping and without another 0.10 cutoff.
  Missing CD8 or CD4 is unavailable. A T cell tied after clipping is reviewed
  across the full identity marker profile plus CD45 before clipping. It keeps
  T-cell identity if CD3e, CD8, CD4 or CD45 is highest; otherwise the ordinary
  chain runs again without T-cell rules. The resolution code records that
  path. Exact preclip ties that remain T cells receive a review record.
- A macrophage is CD68 positive or CD163 positive with CD3e and epithelial
  markers negative. The M1/M2-like subtypes are retired. Each macrophage gets
  a measured CD68/CD163/galectin-3/CD4 profile label if all four values are
  available. Senior's M2Gal3 is the CODEX LGALS3 role; missing values leave
  the profile unavailable rather than supplying a minus sign.
- Compare the resulting hierarchy with a separate preprocessing, UMAP and
  clustering route, as requested for the Barlow-style comparison. Agreement
  is method evidence, not ground truth and not a second production typing
  authority. Do not substitute cluster labels for the hierarchy by default.

These updates are implemented first in the selected pancreas protein chain.
The pLN and spleen pilots still report section-local protein rule candidates
from their recorded panels; neither publishes a final cell class. The owner
explicitly approved the macrophage profile labels on 2026-10-27.

Maintainer clarification, 2026-10-26: tissue-specific typing covers PhenoCycler-only,
Xenium-only and PhenoCycler–Xenium same-slide acquisitions. A protein-only Stage 06
is incomplete. This contract applies to the first pancreas and pancreatic lymph
node release.

## Inputs that Stage 06 must consume

| Input | Producer / authority | Role |
|---|---|---|
| Protein measurements and availability | Stage 03 | Measured protein evidence; absence is not zero |
| Protein normalization | Stage 05 | Donor-consistent marker units and the CD8 shift |
| Xenium cells and RNA QC | Stage 04 (must be built before Stage 06) | RNA evidence on Xenium-only, serial and same-slide sections |
| Tissue, ROI, panel, condition and relationships | Accepted Stage 00 inventory | Route each acquisition without inferring a pair from co-imaging |
| pLN protein phenotype rules | `core/panels/data/tissues/pancreatic_lymph_node/manifest.yml` | Selected full sheet and reviewed pLN adaptation |

The older migration schedule assigned RNA ingest to D8, after D6 typing.
Stage 04 is now a D6 prerequisite: its accepted RNA inputs are needed before
Stage 06 can claim the clarified all-modality typing scope.

Stage 03's measured roster is the 25 pancreas
PhenoCycler images in `cohort-spec-v3.json`. Its accepted receipt does not
cover the 19 dedicated PhenoCycler pLN images. The next measurement increment
must admit those images under the same compartment and exclusion policy before
Stage 06 can call their cells; it must not silently read an older pLN output.

The 2026-10-25 accepted inventory lists **28 donors**: 22 sequential, four
same-slide and two Xenium-only (6512 and 6618). It has 25 PhenoCycler
pancreas images, 19 dedicated PhenoCycler pLN images, 32 Xenium pancreas
regions and 33 Xenium pLN regions. These are acquisitions, not cell counts.
Embedded pLN/lymphoid on a pancreas image must be routed by its recorded
same-block or different-block ROI designation; it is not another acquisition.
The PhenoCycler exclusion for 6457 does not exclude its Xenium data.

The current `io60_v3` cohort spec lists only the 25 PhenoCycler image donors;
it cannot be Stage 04 or Stage 06's sole roster. Stage 00 also records the
Xenium-only 6512/6618 outputs and the four same-slide runs at their
HiPerGator `/orange/maigan` source locations. Those complete outputs are not
on this workstation. Stage 04 must use the accepted acquisition inventory
and resolve those sources on HiPerGator or from a verified transfer. It must
not substitute an older serial Xenium section for a same-slide output.
The current pancreas RNA release worker fixes its roster at 25 sections, and
the pLN v3 RNA worker discovers `/data/xenium/h5ad_out`; that local folder has
neither the two Xenium-only donors nor the four newer same-slide serials.
Reusing either worker's old roster would leave the clarified scope untyped.
The project currently fixes `settings.xenium_root` at local `/data/xenium`;
the relative Xenium bundle paths need a site-specific root on HiPerGator.
The absolute same-slide paths already recorded under `/orange/maigan` must
remain absolute rather than being joined to either root.

`multiplex5000.xenium_sources` reads the accepted Stage 00 acquisition table,
its receipt, and the checksum-bound Xenium records manifest. It uses
`serial__roi` keys, checks complete bundles and their instrument-record hashes,
and reports exact HiPerGator paths for outputs not complete on this host. The
2026-10-25 inventory resolves to 30 ready local region bundles and 35
`remote_only` bundles. `remote_only` means the HiPerGator path is recorded,
not that its large cell files have been verified there; Stage 04 must check
those files on that host before ingest. The tool is a readiness check, not
an RNA producer or a typing result.

`multiplex5000.xenium_region.load_region` is the first Stage 04 reader. It
loads one ready raw bundle's Gene Expression counts and joins its matrix
barcodes to `cells.parquet` by native ID, failing on missing or repeated IDs.
It emits acquisition-qualified `cell_uid` values and keeps cells without a
nucleus in the accounting while marking them ineligible for typing. It makes
no cell-type call or gene-positive cutoff. On the local 6539 pilot, it read
319,005 pancreas cells (41,814 no nucleus) and 222,443 pLN cells (6,344 no
nucleus), both with 5,101 genes. This pilot did not publish Stage 04 or run on
the 35 HPG-only bundles.

The Stage 04 runner declares the location record, accepted Stage 00 evidence,
and each selected region's raw matrix, cell table and instrument records as inputs. It
publishes one region-qualified cell-QC Parquet and feature list per acquisition,
plus a gate figure and receipt; the raw matrices stay at their indexed sources.
Its synthetic same-slide test passes with colliding native IDs, a zero-nucleus
cell, a missing nucleus count, and a damaged pLN region. A project can name
`settings.rna_acquisitions` as an explicit accepted subset for a pilot; the
gate reports its fraction of the full inventory and never claims full coverage.
With no subset it requires all accepted Xenium regions. On the workstation it
can use the accepted local paths and records manifest; on HiPerGator it uses
the complete exact site index. A 65-region HPG run
still needs the accepted Stage 00 evidence staged there and a site file naming
`xenium_source_index`; no full-cohort Stage 04 gate has been accepted.

A read-only HiPerGator check on 2026-10-26 found all seven required files and
matching Stage 00 instrument hashes for all 65 outputs. The 30 outputs that
are complete locally also exist under `/orange/maigan/Xenium_T1D_PancPLN`,
but their exact HiPerGator paths are not all in `MANIFEST.tsv`. The exact
65-region locator is now
`core/multiplex5000/sites/hipergator/xenium_sources_20261026.csv`; Stage 04
must bind it in its own receipt. On HiPerGator all 65 indexed paths passed
the source-file and instrument-hash check. In that tree,
serials 0059865 and 0059907 have pLN outputs named `PanT` rather than the
local `pLN`; the accepted instrument hashes match. Serials 0041116 and
0041117 each have two directory copies whose seven required files compare
byte-identically. The site index should choose one path explicitly, never
pick the first filename match found by a search.

Stage 06 must produce tissue-specific calls, evidence and unavailable reasons
for both modalities. A same-slide section needs both assay results; the
cell-to-cell relationship is established only from the recorded section mode
and registration, then used in Stage 08 integration. Serial sections never
receive a same-cell claim. The damaged 6533 pLN Xenium region may be marked
unavailable with its damage reason and must not hold up other acquisitions.
The current `section_modes.csv` records the four same-slide pancreas serials,
but has no pLN rows for their `pln_2` regions. That gap does not prevent typing
each assay's pLN cells; it prevents a same-cell pLN link until its relationship
is explicitly recorded and checked.

The selected 69-row pLN protein sheet was written for a wider pancreas-plus-pLN
panel. The maintainer approved an explicit pLN adaptation of pancreas-only
negative parent constraints. On the dedicated 55-marker pLN panel, 14 rows
still require unmeasured markers; CD38-dependent B-cell subtype rows remain
unavailable. These protein rules do not type Xenium RNA cells. Stage 06 must
use the tissue-specific RNA producer and map its calls into the release's
class vocabulary. Where the hierarchical protein method writes both
`cell_class` and `refined_class`, only the final step is served.

The pLN and spleen protein pilots report a per-cell **review candidate** from
the sheet's lineage branches: the deepest matched lineage rule only if it is
unique. The panel manifest records that `all` roots are lineages and `any`
roots are orthogonal state annotations. Both retain their rule evidence, but
state matches do not compete for the cell-identity candidate. The pLN manifest
also names the lymph-node lineage roots used for that candidate; pancreas-only
lineages remain in the full sheet as evidence but cannot become pLN identities.
Matched lineage branches remain ambiguous. With no matched lineage,
`no_match` means every candidate lineage root was measured and negative;
`unavailable` means at least one of those roots lacks evidence. A missing
subtype marker alone does not turn a negative root into `unavailable`.
The output records the number of unavailable lineage rules even for a unique
candidate. `candidate_phenotype` is not a
verified final class and cannot be served as one before the processing gate
is accepted.

Stage 06 is the first manual processing review stop. It must show per-tissue,
per-assay and per-acquisition counts, marker or gene support, controls and
unavailable calls. The `NO_VERIFIED_TYPING` latch stays down until the
maintainer accepts processing. No H1–H6 inference runs before that decision.

### Research comparisons and release gates (maintainer clarification, 2026-10-28)

The preceding hold governs **release-backed H1–H6 claims and verified typing**. It
does not bar separate research comparisons on available Xenium data. The HPG
analysis may compute and review exploratory comparative estimates, including
inferential models, before Stage 06 is accepted. Label their inputs and results
as provisional, record the acquisition and donor scope, and do not present a
candidate cell class as verified or activate a release on that basis.

Raw gene-signature comparisons do not depend on Stage 06. A comparison by cell
class needs Stage 06 output, but can read its candidate calls before that gate
is accepted. In the release workflow, Stage 04 acceptance is required to *run*
Stage 06; neither gate blocks a separate analysis of already available source
data. Stages 09 and 10 are not built, so accepting Stage 06 does not start or
unlock a comparative-analysis stage. For the Xenium-only 0041323/0041326 pilot,
the project selects two donors; those two alone cannot support a population-level
disease-effect claim, even if many cells are measured. More donors can enter a
separate, explicitly recorded comparison cohort.

## Spleen in the first portable pilot (maintainer, 2026-10-26)

The first **workflow pilot** also includes executable spleen same-slide typing;
this does not add spleen to the first active IO60 release. The spleen protein
reference is `core/panels/data/tissues/spleen/manifest.yml`. It derives from
the selected full 6551 hierarchy, removes pancreas-only branches and object
labels, retains applicable lymphoid rules, and adds descriptive CD21-positive
B, CD45RO-positive T, CD141-positive cDC and SMA-positive stromal leaves. The
PC67 inventory is the exact, SHA256-checked header from all eight export files.
Of the 58 retained rules, 55 use measured PC67 markers; Neural, Schwann cell
and Macrophage require B3TUBB or Iba1 and stay unavailable. These counts are
rule availability, **not typed-cell counts or validated positivity cutoffs**.
The HPG UniFORM scale table, SHA256 `7927e33b…`, covers only 33 classifier
channels per PC67 section, not the full 67 measured channels. Only five rules
have all inherited markers in that scaled set (Proliferating, PD-L1 expressing,
MHC I high, MHC II high and CD44 high); none establishes a T-cell, B-cell or
FRC lineage. Stage 05 must derive and check the remaining required marker
scales, or mark those lineages unavailable. Raw means in `cells_full` cannot
be substituted for normalized marker evidence. The existing scale table is
comparison evidence, not the completed Stage 05 producer.

The paired protein images and crop-local cell exports are under
`/blue/maigan/smith6jt/SpleenFollicleCounterQP/` on HiPerGator; the eight
Xenium regions are under
`/orange/maigan/20260826__224549__SpleenSH2B3postPheno8sample/`.
The project is an active read-only input for this pilot and is not modified.
Each Xenium region and its corresponding PC67 section must enter the accepted
inventory with a recorded same-slide relationship and coordinate-frame offset.
The two assays receive their own Stage 06 evidence before Stage 08 links cells;
the matching donor code alone never establishes a cell pair. The older spleen
integration plan said RNA cell typing would not run. The maintainer's newer
all-modality typing direction supersedes that scope statement.

The spleen pilot project now routes all eight PC67 full-cell QuPath exports
through Stage 03 before the selected Xenium region enters Stage 04. This Stage
03 adapter preserves `object_id`, uses an acquisition- and export-qualified
`cell_uid`, records all cell means as raw measurements, and excludes zero or
missing nucleus area from protein typing eligibility. It verifies the export
against the accepted Stage 00 digest and labels its gate as raw measurement
evidence. It does not claim the pancreas REDSEA path or normalized protein
values for spleen. Stage 05 needs the eight protein samples to estimate
cross-section shifts; the first bounded RNA pilot uses HDL090's Xenium region.
The other seven Xenium regions remain in the inventory until admitted.

Stage 05 now fits one rigid UniFORM background-mode shift per marker across
the eight Stage 03 PC67 sections. It uses the pinned UniFORM commit through a
site-configured `uniform_checkout` (or a tagged checkout's `external/UniFORM`),
publishes section-by-marker scale and availability tables, and shows before/
after mode spread and individual marker figures. CD8 must have a shift on all
eight sections. A marker with too little measurable signal is reported as
unavailable rather than assigned a synthetic zero histogram. Stage 05 does
not set positivity cutoffs or call cells. Its gate needs processing review
before Stage 06 may use its scales.

**Spleen primary cohort (2026-10-31).** The 12-section `spleen_primary` project replaces
the pilot. It uses the production Protein Measurement (REDSEA → expression) and a
per-panel Normalization route, `stages/panel_normalization.py`, selected by
`settings.normalization: per_panel`:
- UniFORM is fitted on Protein Measurement expression, separately for each panel:
  PC67's eight sections together, PC33's four together, never mixed.
- A marker with fewer than two measurable sections in a panel, or a refused fit, is
  NaN in the aligned cells.
- Its gate is automatic.
- **For this route, Cell Typing does not consume Normalization units.** The spleen
  typing route reads each marker's REDSEA value against that section's own level, so
  Normalization is a prerequisite for typing, not its input.
- The aligned values serve region definitions, classifier rescaling and display.

`multiplex5000.tissues.evaluate.evaluate_rules` is the shared Boolean
hierarchy evaluator for the selected pLN sheet and its spleen adaptation. It
takes explicit per-cell positive, negative or unavailable marker calls; it
does not derive positivity thresholds from intensities or publish a final
cell class. Missing markers in any ancestor make a rule unavailable. Stage 05
and the tissue-specific Stage 06 producers still have to supply and validate
those marker calls before this can type a real cell.

`multiplex5000.tissues.protein_section.evaluate_protein_section` now supplies
that section-local candidate step for both pLN and spleen: measured values,
their Stage 05 shifts, and nucleus admission become tri-state marker calls,
then the tissue's recorded hierarchy is evaluated. It reports missing marker
requirements and makes every rule unavailable for a cell without a nucleus.
The spleen and dedicated pLN Stage 06 pilots use this shared function. The
release-wide pancreas and pLN routes still need the canonical Stage 03 signal
path and all-modality routing; this function is review evidence and does not
publish a final cell class.

The bounded `pln_protein_pilot` project now runs the two recorded dedicated
pLN exports 6380 and 6436 through Stage 03 raw QuPath ingest, Stage 05
UniFORM shifts and Stage 06's selected 69-rule pLN hierarchy. It reads the
portable v3 export index and each site's exact source root. The Stage 03
table retains native `Object ID`, section identity, missing measurements and
nucleus admission; Stage 06 emits per-cell marker and rule states with
unavailable requirements. A synthetic project exercises all three stages and
their task receipts and gates. This project is a review-only code pilot: it
does not replace the approved own-compartment REDSEA path, type Xenium-only
cells, link same-slide cells, serve final classes or activate a release.

The first executable Stage 06 adapter is a bounded spleen same-slide pilot
on HDL090. It applies the maintainer-selected 0.10 background-to-bright-end
candidate independently to each measured protein marker in that physical
section, runs the same hierarchy over those marker calls, and records raw
Xenium identity-gene detection by the existing fixed-depth predicate. The
protein and RNA cells stay in separate tables. The output is review evidence:
there is no `cell_class`, `refined_class`, verified accuracy estimate or
paired-cell result. This adapter does not complete the pancreas/pLN typing
routes or the other seven spleen Xenium regions.

`multiplex5000.xenium_typing_input` now checks a selected Stage 04 region
against its accepted receipt, exact raw bundle digests, output digests and
acquisition-qualified cell IDs before Stage 06 reads it. The spleen adapter
uses that shared check. It can write a disposable raw-count h5ad with recorded
donor and acquisition metadata for the existing pancreas and pLN RNA workers;
that conversion does not select a worker, run its method or publish a class.
The IO60 Stage 06 assay/tissue route still has to call those workers for each
selected acquisition and verify their per-cell outputs.

Both packaged pancreas and pLN RNA workers now accept an exact `cells.parquet`
path in their section functions, with an isolated pLN output root and recorded
ROI. The exact path is checked against the worker h5ad's native cell IDs; a
missing nucleus measurement stops that acquisition's worker call rather than
becoming zero or triggering a donor-level fallback. The Stage 06 route still
has to invoke these functions and report such an acquisition as unavailable
without stopping other acquisitions.
The r4 pancreas uncertainty decisions remain with the frozen r4 release.
They name particular old partitions and cannot be copied onto freshly
clustered acquisitions by donor or section name. A new Stage 06 run uses the
same repaired naming policy but reports that no frozen-partition override was
applied unless the exact partition identity and membership are verified.

## One Stage 06 workflow for both assays (2026-10-28)

`io60_joint_typing_pilot.yaml` is a bounded, portable pancreas+pLN workflow
project. It selects donor 6539 for the packaged mixed-slide protein method and
the two recorded 6539 Xenium regions for the packaged pancreas and pLN RNA
workers. Its Stage 06 route invokes both existing methods in one task, declares
their inputs under one release, and publishes one receipt and gate with separate
protein and RNA checks and figures. Protein and RNA cells keep their own keys;
the route creates no same-cell links before registration and integration. It
does not require typing every acquisition to exercise the workflow.

The project uses the same release for Stage 00, 03, 04, 05 and 06. The site must
name the selected QuPath, geometry, normalization and Xenium inputs. The joint
route reads its own release even if an older site file still names an
`io60_upstream_release` for a separate pilot. From an installed MultiPlex5000 package, use
`mplx5k graph --project io60_joint_typing_pilot --profile local` or select the
HiPerGator profile. The protein hierarchy currently uses the committed Senior
method. The Islet-Explorer session is still refining its CODEX protein rules;
those rules must be reviewed and adopted in this package when that session
finishes. The joint pilot's gate remains a review candidate and does not
activate a release.

The Stage 06 RNA gate checks that every tissue with an intact acquisition that
entered a typing worker has at least one specific lineage call. A processed
tissue with only broad or uncalled cells cannot appear as a successful typing
pilot. Recorded damaged or low-signal acquisitions remain unavailable and do
not block an intact tissue; their counts and reasons stay in the gate report.

RNA clustering also has a recorded numerical setting: `rna_compute_threads: 4`
in each IO60 RNA project, matching the HiPerGator profile's four CPUs per task.
On the same accepted 6539 inputs, two four-thread
wheel runs made identical per-cell calls, while the machine's 64-thread default
made a different but internally repeatable partition. All 28 declared input
hashes matched. Stage 06 now preloads its numerical libraries, limits Numba,
BLAS and OpenMP to the project count while both tissue workers run, and records
library versions and observed thread pools in its receipt. Changing this
setting changes the candidate method and requires a new run; a passing gate
does not make a changed partition biologically correct. Stage 06 refuses to
reuse an existing output whose recorded settings differ from the selected
project; `--force` creates a new candidate and moves the previous output aside.

`io60_xenium_pilot.yaml` now selects the two accepted 6539 Xenium acquisitions
under the existing IO60 candidate release. Its Stage 06 route takes the Stage
04 cells and exact raw source, makes a disposable worker h5ad, and calls the
packaged pancreas or pLN method by the recorded tissue. It verifies that the
per-cell output retains every acquisition-qualified `cell_uid`, counts candidate
and unavailable calls by acquisition, and leaves protein and paired-cell
counts at zero. A recorded damaged pLN section is reported unavailable without
opening its missing raw files or stopping other sections. The synthetic
Stage 04→06 test exercises that routing and receipt/gate path; the heavy
clustering methods are mocked in that test, so a real pilot and scientific
review remain necessary before this route can claim a validated RNA result.
